#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations

import csv, html, json, math, os, re, sys, threading, traceback
from collections import defaultdict
from pathlib import Path
from xml.sax.saxutils import escape as xml_escape
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext

try:
    from PIL import Image, ImageDraw, ImageFont
    PIL_AVAILABLE = True
except Exception:
    Image = ImageDraw = ImageFont = None
    PIL_AVAILABLE = False

VERSION = "5.7-wide-bridge-labels"

def norm_gene(x):
    if x is None: return ""
    s = str(x).strip()
    if not s: return ""
    if "@" in s: s = s.split("@",1)[0]
    return s.lstrip("-").strip()

def norm_interval(x): return str(x or "").strip() or "unknown"

def uniq(xs):
    out, seen = [], set()
    for x in xs:
        k = str(x)
        if k not in seen:
            seen.add(k); out.append(x)
    return out

def interval_bounds(label):
    nums = re.findall(r"\d+(?:\.\d+)?", str(label or ""))
    if len(nums) < 2: return None, None
    try: return float(nums[0]), float(nums[1])
    except Exception: return None, None

def interval_sort_key(label):
    a,b = interval_bounds(label)
    return (10**9,10**9,label) if a is None else (a,b if b is not None else a,label)

def edge_key(a,b): return f"{norm_gene(a).upper()}->{norm_gene(b).upper()}"

def read_tsv(path):
    if not path or not os.path.isfile(path): return []
    with open(path,"r",encoding="utf-8-sig",errors="replace",newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))

def clean_cell_grade(text):
    t = html.unescape(re.sub(r"<[^>]+>"," ",text or ""))
    t = re.sub(r"\s+"," ",t).strip().upper()
    if "NOT EXPRESSED" in t: return "NOT_EXPRESSED"
    if "UNLIKELY" in t: return "UNLIKELY"
    if re.search(r"\bHIGH\b",t): return "HIGH"
    if re.search(r"\bMEDIUM\b",t): return "MEDIUM"
    if re.search(r"\bLOW\b",t): return "LOW"
    return "UNKNOWN"

def parse_cell_validator_html(path):
    if not path or not os.path.isfile(path): return {}
    text = Path(path).read_text(encoding="utf-8",errors="replace")
    row_re = re.compile(r"<tr[^>]*>(?P<row>.*?)</tr>",re.I|re.S)
    td_re = re.compile(r"<td[^>]*>(.*?)</td>",re.I|re.S)
    out={}
    def plain(x):
        x=re.sub(r"<br\s*/?>"," | ",x,flags=re.I)
        x=re.sub(r"<[^>]+>"," ",x)
        return re.sub(r"\s+"," ",html.unescape(x)).strip()
    for rm in row_re.finditer(text):
        row=rm.group("row")
        if "→" not in html.unescape(row): continue
        cells=[plain(x) for x in td_re.findall(row)]
        if len(cells)<10: continue
        a,arrow,b=norm_gene(cells[0]),cells[1],norm_gene(cells[2])
        if not a or not b or "→" not in arrow: continue
        out[edge_key(a,b)] = {
            "gene_a":a,"gene_b":b,
            "expr_a":cells[3] if len(cells)>3 else "",
            "expr_b":cells[4] if len(cells)>4 else "",
            "location_a":cells[5] if len(cells)>5 else "",
            "location_b":cells[6] if len(cells)>6 else "",
            "colocalization":cells[7] if len(cells)>7 else "",
            "string_score":cells[8] if len(cells)>8 else "",
            "grade":clean_cell_grade(cells[9] if len(cells)>9 else ""),
            "evidence":cells[10] if len(cells)>10 else "",
        }
    return out

def parse_phosphosite_tsv(path):
    out={}
    for r in read_tsv(path):
        a,b=norm_gene(r.get("gene_a") or ""),norm_gene(r.get("gene_b") or "")
        if not a or not b:
            e=str(r.get("edge") or "")
            if "->" in e:
                aa,bb=e.split("->",1); a,b=norm_gene(aa),norm_gene(bb)
        if a and b:
            out[edge_key(a,b)]=dict(r)
    return out

def load_temporal_payload(path):
    with open(path,"r",encoding="utf-8-sig") as f: d=json.load(f)
    if not isinstance(d,dict) or not isinstance(d.get("paths"),list):
        raise ValueError("Selected JSON does not contain a Stage 4 temporal 'paths' list.")
    return d

def q_from_target_role(target_role):
    tr=str(target_role or "").lower()
    if "support" in tr: return "supporting"
    if "oppos" in tr: return "opposing"
    if "mixed" in tr or "ambig" in tr: return "mixed"
    return "unresolved"

def _interval_times(intervals):
    if not intervals: return []
    bounds=[]
    for iv in intervals:
        a,b=interval_bounds(iv)
        if a is None: return [f"T{i}" for i in range(len(intervals)+1)]
        bounds.append((a,b))
    vals=[bounds[0][0]]+[b for a,b in bounds]
    def fmt(v):
        if v is None: return "?"
        return f"{int(round(v))}min" if abs(v-round(v))<1e-9 else f"{v:g}min"
    return [fmt(v) for v in vals]

def build_continuous_model(payload,cell_edges,phospho_edges):
    paths=payload.get("paths") or []
    intervals=sorted(uniq(norm_interval(s.get("interval") or s.get("id"))
                          for p in paths for s in (p.get("segments") or [])),key=interval_sort_key)
    iv_index={iv:i for i,iv in enumerate(intervals)}
    boundaries=_interval_times(intervals)
    pool=payload.get("condition_functional_pool") or {}
    excluded={norm_gene(x).upper() for x in (pool.get("excluded_low_activity_genes") or []) if norm_gene(x)}
    nodes,edges,path_map={},{},{}
    def anchor(bi,gene,pid,source=False,target_role="",sites=None):
        k=("anchor",bi,gene)
        if k not in nodes:
            nodes[k]={"key":k,"kind":"anchor","gene":gene,"column":bi*2,"roles":set(),
                      "q_states":set(),"site_tokens":[],"path_ids":set(),
                      "condition_pool_excluded":gene.upper() in excluded}
        n=nodes[k]; n["path_ids"].add(pid)
        if source: n["roles"].add("source_anchor")
        if target_role:
            n["roles"].add("target_anchor"); n["q_states"].add(q_from_target_role(target_role))
        if sites: n["site_tokens"]=uniq(n["site_tokens"]+list(sites))
        return k
    def bridge(ii,iv,gene,pid):
        k=("bridge",ii,gene)
        if k not in nodes:
            nodes[k]={"key":k,"kind":"bridge","gene":gene,"column":ii*2+1,"roles":{"prior_bridge"},
                      "q_states":set(),"site_tokens":[],"path_ids":set(),
                      "condition_pool_excluded":gene.upper() in excluded}
        nodes[k]["path_ids"].add(pid); return k
    def add_edge(a,b,iv,ii,pid,depth):
        ak=(a,b,iv)
        ag,bg=nodes[a]["gene"],nodes[b]["gene"]
        if ak not in edges:
            edges[ak]={"key":ak,"from":a,"to":b,"source_gene":ag,"target_gene":bg,
                       "interval":iv,"interval_index":ii,"path_ids":set(),"count":0,"bridge_depth":depth,
                       "cell":cell_edges.get(edge_key(ag,bg),{}),
                       "phospho":phospho_edges.get(edge_key(ag,bg),{})}
        edges[ak]["path_ids"].add(pid); edges[ak]["count"]+=1

    for p_idx,p in enumerate(paths,1):
        pid=str(p.get("path_id") or p.get("id") or p_idx)
        path_map[pid]={"path_id":pid,"genes":[norm_gene(x) for x in (p.get("genes") or []) if norm_gene(x)],
                       "node_keys":[],"edge_keys":[]}
        for seg in p.get("segments") or []:
            iv=norm_interval(seg.get("interval") or seg.get("id")); ii=iv_index.get(iv,0)
            src,tgt=norm_gene(seg.get("source_anchor")),norm_gene(seg.get("target_anchor"))
            bridges=[norm_gene(x) for x in (seg.get("bridge_nodes") or []) if norm_gene(x)]
            if not src or not tgt: continue
            sk=anchor(ii,src,pid,source=True)
            tk=anchor(ii+1,tgt,pid,target_role=str(seg.get("target_role") or "unresolved"),
                      sites=[str(x) for x in (seg.get("target_anchor_tokens") or []) if str(x).strip()])
            chain=[sk]+[bridge(ii,iv,b,pid) for b in bridges]+[tk]
            path_map[pid]["node_keys"].extend(chain)
            for j in range(len(chain)-1):
                add_edge(chain[j],chain[j+1],iv,ii,pid,len(bridges))
                path_map[pid]["edge_keys"].append((chain[j],chain[j+1],iv))
    for n in nodes.values():
        qs=n["q_states"]
        if n["kind"]=="bridge": n["q_state"]="na"
        elif len(qs)>1: n["q_state"]="mixed"
        elif len(qs)==1: n["q_state"]=next(iter(qs))
        else: n["q_state"]="source"
        n["roles"]=sorted(n["roles"]); n["q_states"]=sorted(n["q_states"])
    columns=defaultdict(list)
    for k,n in nodes.items(): columns[n["column"]].append(k)
    order={}
    for c in sorted(columns):
        columns[c]=sorted(columns[c],key=lambda k:nodes[k]["gene"])
        for i,k in enumerate(columns[c]): order[k]=i
    L,R=defaultdict(list),defaultdict(list)
    for e in edges.values():
        a,b=e["from"],e["to"]
        if nodes[a]["column"]<=nodes[b]["column"]: R[a].append(b); L[b].append(a)
        else: R[b].append(a); L[a].append(b)
    def reorder(c,use_left):
        arr=columns[c]; scored=[]
        for k in arr:
            nb=L[k] if use_left else R[k]
            bary=(sum(order.get(x,0) for x in nb)/len(nb)) if nb else order.get(k,0)
            scored.append((bary,nodes[k]["gene"],k))
        scored.sort(key=lambda x:(x[0],x[1])); columns[c]=[x[2] for x in scored]
        for i,k in enumerate(columns[c]): order[k]=i
    cols=sorted(columns)
    for _ in range(8):
        for c in cols[1:]: reorder(c,True)
        for c in reversed(cols[:-1]): reorder(c,False)
    for c in cols:
        for i,k in enumerate(columns[c]): nodes[k]["lane_index"]=i
    matched_cell_edges=sum(1 for e in edges.values() if e.get("cell"))
    matched_phospho_edges=sum(1 for e in edges.values() if e.get("phospho"))
    return {"version":VERSION,"intervals":intervals,"boundary_labels":boundaries,"nodes":list(nodes.values()),
            "edges":list(edges.values()),"paths":path_map,"pool":pool,
            "bridge_depth_D":payload.get("bridge_depth_D"),"max_edge_count":max([e["count"] for e in edges.values()] or [1]),
            "validation_summary":{"cell_records_loaded":len(cell_edges),"phospho_records_loaded":len(phospho_edges),
                                  "cell_edges_matched":matched_cell_edges,"phospho_edges_matched":matched_phospho_edges}}

class ContinuousCanvas(tk.Canvas):
    COLUMN_GAP=182; TOP=95; ROW_GAP=44; NODE_H=26; NODE_PAD_X=12
    def __init__(self,master,detail_callback,**kw):
        super().__init__(master,background="white",highlightthickness=0,**kw)
        self.model=None; self.detail_callback=detail_callback
        self.scale_factor=1.0; self.pan_x=40; self.pan_y=25; self.drag_start=None
        self.node_item_to_record={}; self.edge_item_to_record={}
        self.selected_path=""; self.search_gene=""
        self.enabled_q={"supporting","opposing","unresolved","mixed","source","na"}

        # Floating legend positions are stored as fractions of the current canvas
        # so their placement is preserved when resizing and can be reused in PNG/SVG export.
        self.legend_pos={"top":[0.06,0.035], "bottom":[0.06,0.82]}
        self.legend_drag=None

        self.bind("<ButtonPress-1>",self._press); self.bind("<B1-Motion>",self._drag); self.bind("<ButtonRelease-1>",self._release)
        self.bind("<MouseWheel>",self._wheel); self.bind("<Double-Button-1>",lambda e:self.fit_overview())
    def set_model(self,m): self.model=m; self.readable_view()
    def readable_view(self): self.scale_factor=1.0; self.pan_x=45; self.pan_y=25; self.redraw()
    def zoom_in(self): self._zoom_at(self.winfo_width()/2,self.winfo_height()/2,1.18)
    def zoom_out(self): self._zoom_at(self.winfo_width()/2,self.winfo_height()/2,0.85)
    def spread_vertical(self,f): self.ROW_GAP=int(max(30,min(90,self.ROW_GAP*f))); self.redraw()
    def spread_horizontal(self,f): self.COLUMN_GAP=int(max(120,min(320,self.COLUMN_GAP*f))); self.redraw()
    def fit_overview(self):
        if not self.model:return
        self.after_idle(self._fit)
    def _fit(self):
        max_col=max([n["column"] for n in self.model["nodes"]] or [0]); max_rows=self._max_col_count()
        ww=max(900,max_col*self.COLUMN_GAP+320); wh=max(600,self.TOP+max_rows*self.ROW_GAP+120)
        self.scale_factor=min(max(100,self.winfo_width()-40)/ww,max(100,self.winfo_height()-40)/wh,1.0)
        self.pan_x=25; self.pan_y=20; self.redraw()
    def _world_to_screen(self,x,y): return x*self.scale_factor+self.pan_x,y*self.scale_factor+self.pan_y
    def _screen_to_world(self,x,y): return (x-self.pan_x)/self.scale_factor,(y-self.pan_y)/self.scale_factor
    def _max_col_count(self):
        d=defaultdict(int)
        for n in self.model["nodes"]: d[n["column"]]+=1
        return max(d.values() or [1])
    def _node_world_pos(self,n):
        x=120+n["column"]*self.COLUMN_GAP
        coln=[m for m in self.model["nodes"] if m["column"]==n["column"]]
        y0=self.TOP+(self._max_col_count()-len(coln))*self.ROW_GAP/2
        return x,y0+n.get("lane_index",0)*self.ROW_GAP
    def _bridge_polygon_points(self, x, y, width, height):
        # Horizontally stretched diamond / hexagon so bridge labels remain readable in print.
        tip = max(10, min(width * 0.18, 18 * self.scale_factor))
        return [
            x - width/2 + tip, y - height/2,
            x + width/2 - tip, y - height/2,
            x + width/2,       y,
            x + width/2 - tip, y + height/2,
            x - width/2 + tip, y + height/2,
            x - width/2,       y,
        ]

    def _visible_nodes(self):
        pa=None
        if self.selected_path:
            p=self.model["paths"].get(self.selected_path); pa=set(p.get("node_keys",[])) if p else set()
        out=set()
        for n in self.model["nodes"]:
            k=tuple(n["key"])
            if n["q_state"] not in self.enabled_q: continue
            if pa is not None and k not in pa: continue
            out.add(k)
        return out
    def _node_style(self,n):
        if n["kind"]=="bridge": return "#ffffff","#2563eb","wide_bridge"
        return {
            "supporting":("#22c55e","#166534","round"),
            "opposing":("#ef4444","#991b1b","round"),
            "mixed":("#8b5cf6","#5b21b6","round"),
            "unresolved":("#f59e0b","#92400e","round"),
            "source":("#1d4ed8","#1e3a8a","round"),
        }.get(n["q_state"],("#e2e8f0","#475569","round"))
    def _edge_style(self,e):
        ph=str((e.get("phospho") or {}).get("grade","")).upper()
        cg=str((e.get("cell") or {}).get("grade","UNKNOWN")).upper()
        nl={tuple(n["key"]):n for n in self.model["nodes"]}
        q=nl.get(tuple(e["to"]),{}).get("q_state","unresolved")
        if ph.startswith(("A","B")): color="#166534"
        elif ph.startswith("C"): color="#d97706"
        elif cg in ("UNLIKELY","NOT_EXPRESSED"): color="#dc2626"
        else: color={"supporting":"#15803d","opposing":"#dc2626","mixed":"#7c3aed","unresolved":"#d97706"}.get(q,"#475569")
        m=max(1,self.model.get("max_edge_count",1)); width=1+6*((e.get("count",1)/m)**0.60)
        return color,width
    def redraw(self):
        self.delete("all"); self.node_item_to_record.clear(); self.edge_item_to_record.clear()
        if not self.model:
            self.create_text(30,30,anchor="nw",text="Select a Stage 4.1 temporal JSON.",fill="#64748b"); return
        vis=self._visible_nodes(); max_col=max([n["column"] for n in self.model["nodes"]] or [0])
        bounds=self.model.get("boundary_labels") or []; ivs=self.model.get("intervals") or []
        maxwy=self.TOP+self._max_col_count()*self.ROW_GAP+120
        for c in range(max_col+1):
            wx=120+c*self.COLUMN_GAP; x,_=self._world_to_screen(wx,0); _,y0=self._world_to_screen(0,10); _,y1=self._world_to_screen(0,maxwy)
            fill="#eff6ff" if c%2==0 else "#fffaf0"; half=self.COLUMN_GAP*0.46*self.scale_factor
            self.create_rectangle(x-half,y0,x+half,y1,fill=fill,outline="#e2e8f0")
            if c%2==0:
                bi=c//2; label=bounds[bi] if bi<len(bounds) else f"T{bi}"
            else:
                ii=c//2; label=f"Bridge\n{ivs[ii] if ii<len(ivs) else ''}"
            self.create_text(x,y0+20,text=label,font=("Segoe UI",10,"bold"),fill="#0f172a",justify="center")
        lookup={tuple(n["key"]):n for n in self.model["nodes"]}
        for e in self.model["edges"]:
            a,b=tuple(e["from"]),tuple(e["to"])
            if a not in vis or b not in vis: continue
            if self.selected_path and self.selected_path not in e["path_ids"]: continue
            x1,y1=self._world_to_screen(*self._node_world_pos(lookup[a])); x2,y2=self._world_to_screen(*self._node_world_pos(lookup[b]))
            col,w=self._edge_style(e); mx=(x1+x2)/2; my=(y1+y2)/2+(((hash((e["source_gene"],e["target_gene"],e["interval"]))%7)-3)*2.0)
            item=self.create_line(x1,y1,mx,my,x2,y2,fill=col,width=max(1,w*self.scale_factor),arrow="last",arrowshape=(9,11,4),smooth=True,splinesteps=16)
            self.edge_item_to_record[item]=e
        for n in self.model["nodes"]:
            k=tuple(n["key"])
            if k not in vis: continue
            x,y=self._world_to_screen(*self._node_world_pos(n)); fill,border,shape=self._node_style(n)
            hit=bool(self.search_gene and self.search_gene in n["gene"].upper())
            tw=max(44,len(n["gene"])*8+24)*self.scale_factor
            h=self.NODE_H*self.scale_factor
            if shape=="wide_bridge":
                b_w=max(66, len(n["gene"])*8+28)*self.scale_factor
                pts=self._bridge_polygon_points(x, y, b_w, h)
                item=self.create_polygon(pts,fill=fill,outline="#111827" if hit else border,width=4 if hit else 2)
            else:
                item=self.create_oval(x-tw/2,y-h/2,x+tw/2,y+h/2,fill=fill,outline="#111827" if hit else border,width=4 if hit else 2)
            self.node_item_to_record[item]=n
            fs=max(8,min(13,int(10.5*self.scale_factor))); tc="white" if n["q_state"]=="source" else "#111827"
            t=self.create_text(x,y,text=n["gene"],font=("Segoe UI",fs,"bold"),fill=tc); self.node_item_to_record[t]=n

        self._draw_overlay_legends()
        self.configure(scrollregion=self.bbox("all") or (0,0,1000,800))
    def reset_legend_positions(self):
        self.legend_pos={"top":[0.06,0.035], "bottom":[0.06,0.82]}
        self.redraw()

    def _press(self,e):
        cur=self.find_withtag("current")
        if cur:
            i=cur[0]
            tags=self.gettags(i)
            if "legend_top" in tags:
                self.legend_drag=("top",e.x,e.y)
                return
            if "legend_bottom" in tags:
                self.legend_drag=("bottom",e.x,e.y)
                return
            if i in self.node_item_to_record:
                self.detail_callback("node",self.node_item_to_record[i]); return
            if i in self.edge_item_to_record:
                self.detail_callback("edge",self.edge_item_to_record[i]); return
        self.drag_start=(e.x,e.y)

    def _drag(self,e):
        if self.legend_drag:
            which,lastx,lasty=self.legend_drag
            w=max(1,self.winfo_width()); h=max(1,self.winfo_height())
            dx=(e.x-lastx)/w; dy=(e.y-lasty)/h
            fx,fy=self.legend_pos[which]
            # Keep the legend mostly inside the canvas.
            fx=max(0.0,min(0.80,fx+dx))
            fy=max(0.0,min(0.90,fy+dy))
            self.legend_pos[which]=[fx,fy]
            self.legend_drag=(which,e.x,e.y)
            self.redraw()
            return
        if not self.drag_start:return
        self.pan_x+=e.x-self.drag_start[0]; self.pan_y+=e.y-self.drag_start[1]
        self.drag_start=(e.x,e.y); self.redraw()

    def _release(self,e):
        self.drag_start=None
        self.legend_drag=None

    def _wheel(self,e):
        ctrl=bool(e.state & 0x0004); shift=bool(e.state & 0x0001)
        if ctrl:self._zoom_at(e.x,e.y,1.12 if e.delta>0 else 0.89)
        elif shift:self.pan_x+=80 if e.delta>0 else -80; self.redraw()
        else:self.pan_y+=90 if e.delta>0 else -90; self.redraw()
    def _zoom_at(self,sx,sy,f):
        wx,wy=self._screen_to_world(sx,sy); new=max(0.18,min(3.0,self.scale_factor*f))
        self.scale_factor=new; self.pan_x=sx-wx*new; self.pan_y=sy-wy*new; self.redraw()

    def _draw_overlay_legends(self):
        """Draggable floating legends. Drag either legend by any point inside its box."""
        w = max(1000, self.winfo_width())
        h = max(700, self.winfo_height())
        vs = self.model.get("validation_summary") or {}
        cell_loaded = int(vs.get("cell_records_loaded",0))
        ph_loaded = int(vs.get("phospho_records_loaded",0))
        cell_matched = int(vs.get("cell_edges_matched",0))
        ph_matched = int(vs.get("phospho_edges_matched",0))

        # ---------- TOP NODE / Q LEGEND ----------
        fx,fy=self.legend_pos["top"]
        x0=max(6,fx*w); y0=max(6,fy*h)
        box_w=min(1180,max(760,w*0.76))
        box_h=100
        x1=min(w-6,x0+box_w); y1=min(h-6,y0+box_h)
        tag=("legend_top","floating_legend")
        self.create_rectangle(x0,y0,x1,y1,fill="#ffffff",outline="#64748b",width=2,tags=tag)
        self.create_text(x0+14,y0+10,anchor="nw",text="NODE / Q-STATE LEGEND   [drag to move]",
                         font=("Segoe UI",11,"bold"),fill="#0f172a",tags=tag)
        items=[
            ("circle","#1d4ed8","#1e3a8a","Initial / boundary source"),
            ("circle","#22c55e","#166534","Supporting anchor"),
            ("circle","#f59e0b","#92400e","Unresolved anchor"),
            ("circle","#8b5cf6","#5b21b6","Mixed anchor"),
            ("circle","#ef4444","#991b1b","Opposing anchor"),
            ("wide_bridge","#ffffff","#2563eb","Prior-knowledge bridge"),
        ]
        usable=max(300,x1-x0-28); step=usable/3.0
        for i,(shape,fill,border,label) in enumerate(items):
            row=i//3; col=i%3
            bx=x0+14+col*step; by=y0+43+row*29
            cx=bx+10; cy=by+7
            if shape=="wide_bridge":
                ww,hh=34,18; tip=8
                self.create_polygon([cx-ww/2+tip,cy-hh/2,cx+ww/2-tip,cy-hh/2,cx+ww/2,cy,cx+ww/2-tip,cy+hh/2,cx-ww/2+tip,cy+hh/2,cx-ww/2,cy],
                                    fill=fill,outline=border,width=2,tags=tag)
            else:
                self.create_oval(cx-9,cy-9,cx+9,cy+9,fill=fill,outline=border,width=2,tags=tag)
            self.create_text(bx+28,cy,anchor="w",text=label,
                             font=("Segoe UI",10,"bold"),fill="#111827",tags=tag)

        # ---------- BOTTOM EDGE / VALIDATION LEGEND ----------
        fx,fy=self.legend_pos["bottom"]
        x0b=max(6,fx*w); y0b=max(6,fy*h)
        box_w=min(1330,max(860,w*0.84))
        box_h=128
        x1b=min(w-6,x0b+box_w); y1b=min(h-6,y0b+box_h)
        tagb=("legend_bottom","floating_legend")
        self.create_rectangle(x0b,y0b,x1b,y1b,fill="#ffffff",outline="#64748b",width=2,tags=tagb)
        self.create_text(x0b+14,y0b+10,anchor="nw",text="EDGE / VALIDATION LEGEND   [drag to move]",
                         font=("Segoe UI",11,"bold"),fill="#0f172a",tags=tagb)
        specs=[
            ("#15803d",4,"Green: supporting or phosphosite Grade A/B"),
            ("#d97706",4,"Orange: unresolved or phosphosite Grade C"),
            ("#dc2626",4,"Red: opposing or cell-type UNLIKELY / NOT_EXPRESSED"),
            ("#475569",2,"Gray-blue: general functional association"),
        ]
        usable=max(400,x1b-x0b-28); step=usable/2.0
        for i,(color,width,label) in enumerate(specs):
            row=i//2; col=i%2
            bx=x0b+14+col*step; by=y0b+42+row*27
            self.create_line(bx,by,bx+44,by,fill=color,width=width,
                             arrow="last",arrowshape=(9,11,4),tags=tagb)
            self.create_text(bx+56,by,anchor="w",text=label,
                             font=("Segoe UI",10),fill="#111827",tags=tagb)

        cell_txt=(f"Cell-type validation: LOADED ({cell_loaded} records; {cell_matched} matched graph edges)"
                  if cell_loaded else "Cell-type validation: NOT LOADED")
        ph_txt=(f"Phosphosite validation: LOADED ({ph_loaded} records; {ph_matched} matched graph edges)"
                if ph_loaded else "Phosphosite validation: NOT LOADED")
        self.create_text(x0b+14,y1b-31,anchor="sw",text=cell_txt,
                         font=("Segoe UI",9,"bold"),fill="#334155",tags=tagb)
        self.create_text(x0b+14,y1b-10,anchor="sw",text=ph_txt,
                         font=("Segoe UI",9,"bold"),fill="#334155",tags=tagb)
        self.create_text(x1b-14,y1b-10,anchor="se",
                         text="Edge thickness ∝ repeated route occurrences",
                         font=("Segoe UI",9),fill="#475569",tags=tagb)


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("PPI Temporal Graphics — Stage 5.7 Wide-Bridge Validation-Aware Viewer")
        self.geometry("1720x980")
        self.minsize(1180,720)
        self.fullscreen_on = False
        self.left_panel_visible = True
        self.right_panel_visible = True
        self.temporal_var=tk.StringVar(); self.cell_var=tk.StringVar(); self.phospho_var=tk.StringVar()
        self.status_var=tk.StringVar(value="Ready."); self.path_var=tk.StringVar(value="All paths"); self.search_var=tk.StringVar()
        self.model=None; self.q_vars={}
        self._build()
        try:
            self.state("zoomed")  # Windows maximize
        except Exception:
            pass
        self.bind("<F11>", lambda e: self.toggle_fullscreen())
        self.bind("<Escape>", lambda e: self.exit_fullscreen())
    def _build(self):
        top=ttk.Frame(self,padding=8); top.pack(fill="x")
        ttk.Label(top,text="Stage 5.2 — Continuous Temporal PPI Network + Validation Overlay",font=("Segoe UI",16,"bold")).grid(row=0,column=0,columnspan=12,sticky="w")
        ttk.Label(top,text="Shared temporal-boundary anchors • compact left-to-right flow • repeated-edge weighting • in-canvas legends",foreground="#475569").grid(row=1,column=0,columnspan=12,sticky="w",pady=(0,6))
        ttk.Label(top,text="Temporal JSON").grid(row=2,column=0,sticky="w")
        ttk.Entry(top,textvariable=self.temporal_var,width=75).grid(row=2,column=1,columnspan=6,sticky="ew",padx=5)
        ttk.Button(top,text="Browse...",command=self.choose_temporal).grid(row=2,column=7)
        ttk.Button(top,text="Load / Draw",command=self.load_draw).grid(row=2,column=8,padx=5)
        ttk.Button(top,text="Save SVG",command=self.save_svg).grid(row=2,column=9,padx=(14,3))
        ttk.Button(top,text="Export CSV",command=self.export_csv).grid(row=2,column=10,padx=3)
        ttk.Button(top,text="Open Folder",command=self.open_folder).grid(row=2,column=11,padx=3)
        ttk.Label(top,text="Cell validator HTML").grid(row=3,column=0,sticky="w")
        ttk.Entry(top,textvariable=self.cell_var).grid(row=3,column=1,columnspan=3,sticky="ew",padx=5)
        ttk.Button(top,text="Browse...",command=self.choose_cell).grid(row=3,column=4)
        ttk.Label(top,text="Phosphosite TSV").grid(row=3,column=5,sticky="e")
        ttk.Entry(top,textvariable=self.phospho_var).grid(row=3,column=6,columnspan=4,sticky="ew",padx=5)
        ttk.Button(top,text="Browse...",command=self.choose_phospho).grid(row=3,column=10)
        ctl=ttk.Frame(top); ctl.grid(row=4,column=0,columnspan=12,sticky="w",pady=(7,0))
        for txt,cmd in [
                        ("Readable",lambda:self.canvas.readable_view()),
                        ("Overview",lambda:self.canvas.fit_overview()),
                        ("Zoom +",lambda:self.canvas.zoom_in()),
                        ("Zoom −",lambda:self.canvas.zoom_out()),
                        ("More vertical space",lambda:self.canvas.spread_vertical(1.18)),
                        ("Less vertical space",lambda:self.canvas.spread_vertical(0.85)),
                        ("More horizontal space",lambda:self.canvas.spread_horizontal(1.15)),
                        ("Less horizontal space",lambda:self.canvas.spread_horizontal(0.87)),
                        ("Toggle Fullscreen (F11)",self.toggle_fullscreen),
                        ("Canvas Focus",self.toggle_canvas_focus),
                        ("Reset Legends",lambda:self.canvas.reset_legend_positions()),
                        ]:
            ttk.Button(ctl,text=txt,command=cmd).pack(side="left",padx=3)

        export_ctl=ttk.Frame(top); export_ctl.grid(row=5,column=0,columnspan=12,sticky="w",pady=(4,0))
        ttk.Label(export_ctl, text="Export / restored utilities:", foreground="#475569").pack(side="left", padx=(0,6))
        for txt,cmd in [
                        ("Save PNG",self.save_png),
                        ("Save Hi-Res PNG",self.save_hires_png),
                        ("Save SVG",self.save_svg),
                        ("Export CSV",self.export_csv),
                        ("Open Folder",self.open_folder),
                        ]:
            ttk.Button(export_ctl,text=txt,command=cmd).pack(side="left",padx=3)
        for c in range(12): top.columnconfigure(c,weight=1 if c in (1,2,3,6) else 0)
        self.main=ttk.Panedwindow(self,orient="horizontal"); self.main.pack(fill="both",expand=True,padx=8,pady=(0,8))
        self.left_frame,center,self.right_frame=ttk.Frame(self.main,width=230),ttk.Frame(self.main),ttk.Frame(self.main,width=345)
        self.main.add(self.left_frame,weight=0); self.main.add(center,weight=1); self.main.add(self.right_frame,weight=0)
        left = self.left_frame
        right = self.right_frame
        qf=ttk.LabelFrame(left,text="Q state",padding=6); qf.pack(fill="x",padx=4,pady=4)
        for q,label in [("supporting","Supporting anchor"),("unresolved","Unresolved anchor"),("mixed","Mixed anchor"),("opposing","Opposing anchor"),("source","Initial source"),("na","Prior bridge")]:
            v=tk.BooleanVar(value=True); self.q_vars[q]=v; ttk.Checkbutton(qf,text=label,variable=v,command=self.apply_filters).pack(anchor="w")
        pf=ttk.LabelFrame(left,text="Path / gene",padding=6); pf.pack(fill="x",padx=4,pady=4)
        self.path_combo=ttk.Combobox(pf,textvariable=self.path_var,state="readonly",values=["All paths"]); self.path_combo.pack(fill="x")
        self.path_combo.bind("<<ComboboxSelected>>",lambda e:self.apply_filters())
        ttk.Entry(pf,textvariable=self.search_var).pack(fill="x",pady=(7,0)); self.search_var.trace_add("write",lambda *_:self.apply_filters())
        ttk.Label(pf,text="Gene search highlights; it does not hide the route.",foreground="#64748b",wraplength=200).pack(anchor="w")
        leg=ttk.LabelFrame(left,text="Legend",padding=6); leg.pack(fill="x",padx=4,pady=4)
        for txt,col in [("Initial source","#1d4ed8"),("Supporting","#22c55e"),("Unresolved","#f59e0b"),("Mixed","#8b5cf6"),("Opposing","#ef4444"),("Prior bridge","#ffffff")]:
            r=ttk.Frame(leg); r.pack(fill="x",pady=2); sw=tk.Canvas(r,width=17,height=17,highlightthickness=0); sw.create_oval(2,2,15,15,fill=col,outline="#2563eb" if txt=="Prior bridge" else "#475569"); sw.pack(side="left"); ttk.Label(r,text=txt).pack(side="left",padx=5)
        ttk.Label(leg,text="Edge thickness = repeated route occurrence count.",foreground="#64748b",wraplength=200).pack(anchor="w",pady=(5,0))
        cf=ttk.Frame(center); cf.pack(fill="both",expand=True)
        xb,yb=ttk.Scrollbar(cf,orient="horizontal"),ttk.Scrollbar(cf,orient="vertical")
        self.canvas=ContinuousCanvas(cf,self.show_detail,xscrollcommand=xb.set,yscrollcommand=yb.set)
        xb.configure(command=self.canvas.xview); yb.configure(command=self.canvas.yview)
        self.canvas.grid(row=0,column=0,sticky="nsew"); yb.grid(row=0,column=1,sticky="ns"); xb.grid(row=1,column=0,sticky="ew")
        cf.rowconfigure(0,weight=1); cf.columnconfigure(0,weight=1)
        ttk.Label(center,textvariable=self.status_var,anchor="w").pack(fill="x",pady=(3,0))
        rf=ttk.LabelFrame(right,text="Selected node / edge",padding=6); rf.pack(fill="both",expand=True,padx=4,pady=4)
        self.detail=scrolledtext.ScrolledText(rf,width=43,font=("Consolas",9),state="disabled"); self.detail.pack(fill="both",expand=True)
        pf2=ttk.LabelFrame(right,text="Condition-specific functional pool",padding=6); pf2.pack(fill="x",padx=4,pady=4)
        self.pool_text=tk.Text(pf2,height=7,width=42,font=("Consolas",8),state="disabled"); self.pool_text.pack(fill="x")
        ttk.Label(self,text="Wheel = vertical pan | Shift+Wheel = horizontal pan | Ctrl+Wheel = zoom | drag blank space = pan | double-click = overview | F11 = fullscreen",foreground="#64748b",anchor="w",padding=(8,0,8,7)).pack(fill="x")

    def _default_output_path(self, ext: str, suffix: str = "_stage5"):
        src = self.temporal_var.get().strip()
        if src and os.path.isfile(src):
            p = Path(src)
            return str(p.with_name(p.stem + suffix + ext))
        return str(Path.cwd() / f"temporal_viewer{suffix}{ext}")

    def toggle_fullscreen(self):
        self.fullscreen_on = not self.fullscreen_on
        try:
            self.attributes("-fullscreen", self.fullscreen_on)
        except Exception:
            try:
                self.state("zoomed" if self.fullscreen_on else "normal")
            except Exception:
                pass

    def exit_fullscreen(self):
        self.fullscreen_on = False
        try:
            self.attributes("-fullscreen", False)
        except Exception:
            pass

    def toggle_canvas_focus(self):
        # Hide or restore side panels to maximize drawing area.
        try:
            panes = list(self.main.panes())
            if self.left_panel_visible or self.right_panel_visible:
                if self.left_panel_visible and str(self.left_frame) in panes:
                    self.main.forget(self.left_frame)
                if self.right_panel_visible and str(self.right_frame) in panes:
                    self.main.forget(self.right_frame)
                self.left_panel_visible = False
                self.right_panel_visible = False
                self.status_var.set("Canvas focus mode enabled. Click 'Canvas Focus' again to restore side panels.")
            else:
                self.main.insert(0, self.left_frame, weight=0)
                self.main.add(self.right_frame, weight=0)
                self.left_panel_visible = True
                self.right_panel_visible = True
                self.status_var.set("Side panels restored.")
        except Exception as e:
            messagebox.showerror("Canvas Focus Error", str(e))

    def choose_temporal(self):
        p=filedialog.askopenfilename(title="Select Stage 4.1 temporal JSON",filetypes=[("JSON","*.json"),("All files","*.*")])
        if p:self.temporal_var.set(p)
    def choose_cell(self):
        p=filedialog.askopenfilename(title="Select cell-type validator HTML",filetypes=[("HTML","*.html;*.htm"),("All files","*.*")])
        if p:self.cell_var.set(p)
    def choose_phospho(self):
        p=filedialog.askopenfilename(title="Select phosphosite validation TSV",filetypes=[("TSV","*.tsv"),("All files","*.*")])
        if p:self.phospho_var.set(p)
    def load_draw(self):
        p=self.temporal_var.get().strip()
        if not os.path.isfile(p): messagebox.showerror("Input Error","Please select a valid Stage 4.1 temporal JSON."); return
        self.status_var.set("Loading..."); threading.Thread(target=self._worker,daemon=True).start()
    def _worker(self):
        try:
            pl=load_temporal_payload(self.temporal_var.get().strip())
            ce=parse_cell_validator_html(self.cell_var.get().strip()) if self.cell_var.get().strip() else {}
            pe=parse_phosphosite_tsv(self.phospho_var.get().strip()) if self.phospho_var.get().strip() else {}
            m=build_continuous_model(pl,ce,pe)
            m["validation_summary"]["cell_file_selected"]=bool(self.cell_var.get().strip())
            m["validation_summary"]["phospho_file_selected"]=bool(self.phospho_var.get().strip())
            self.after(0,lambda:self._finish(m))
        except Exception as e:
            tb=traceback.format_exc(); self.after(0,lambda:messagebox.showerror("Stage 5.2 Error",f"{e}\n\n{tb}")); self.after(0,lambda:self.status_var.set("Error."))
    def _finish(self,m):
        self.model=m
        self.path_combo["values"]=["All paths"]+[f"Path {pid}: {' → '.join(r['genes'])}" for pid,r in m["paths"].items()]
        self.path_var.set("All paths"); self.canvas.set_model(m); self._show_pool()
        ac=sum(1 for n in m["nodes"] if n["kind"]=="anchor"); bc=sum(1 for n in m["nodes"] if n["kind"]=="bridge")
        vs=m.get("validation_summary") or {}
        self.status_var.set(f"Loaded: {len(m['paths'])} paths | {ac} anchors | {bc} bridges | {len(m['edges'])} edges | "
                            f"Cell validation matched={vs.get('cell_edges_matched',0)} | Phosphosite validation matched={vs.get('phospho_edges_matched',0)} | "
                            f"intervals: {', '.join(m['intervals'])}.")
    def apply_filters(self):
        if not self.model:return
        self.canvas.enabled_q={q for q,v in self.q_vars.items() if v.get()}
        t=self.path_var.get(); self.canvas.selected_path=""
        if t.startswith("Path "): self.canvas.selected_path=t.split(":",1)[0].replace("Path ","").strip()
        self.canvas.search_gene=self.search_var.get().strip().upper(); self.canvas.redraw()
    def _show_pool(self):
        p=self.model.get("pool") or {}; txt=f"Basis: {p.get('basis','not specified')}\nExcluded low-activity genes: {len(p.get('excluded_low_activity_genes') or [])}\n\n{', '.join(p.get('excluded_low_activity_genes') or []) or 'None'}"
        self.pool_text.configure(state="normal"); self.pool_text.delete("1.0","end"); self.pool_text.insert("1.0",txt); self.pool_text.configure(state="disabled")
    def show_detail(self,kind,r):
        if kind=="node":
            lines=[f"NODE: {r['gene']}",f"Visual kind: {r['kind']}",f"Roles: {', '.join(r.get('roles') or []) or '—'}",f"Q state: {r.get('q_state','—')}",f"Observed sites: {', '.join(r.get('site_tokens') or []) or '—'}",f"Condition-pool excluded: {r.get('condition_pool_excluded')}",f"Path count represented: {len(r.get('path_ids') or [])}","","Paths:",", ".join(sorted(r.get("path_ids") or []))]
        else:
            c=r.get("cell") or {}; p=r.get("phospho") or {}
            lines=[f"EDGE: {r.get('source_gene')} -> {r.get('target_gene')}",f"Interval: {r.get('interval')}",f"Repeated path occurrences: {r.get('count',1)}",f"Bridge depth: {r.get('bridge_depth','—')}","","CELL-TYPE VALIDATION",f"Grade: {c.get('grade','Unknown')}",f"Expression A: {c.get('expr_a','—')}",f"Expression B: {c.get('expr_b','—')}",f"Co-localization: {c.get('colocalization','—')}",f"STRING: {c.get('string_score','—')}","","PHOSPHOSITE EVIDENCE",f"Grade: {p.get('grade','Unknown')}",f"Confidence: {p.get('confidence_score','—')}",f"Matched sites: {p.get('matched_sites','—')}",f"SIGNOR residues: {p.get('signor_residues','—')}",f"SIGNOR effects: {p.get('signor_effects','—')}",f"OmniPath: {p.get('omnipath_enzsub_sources') or p.get('omnipath_interact_sources') or '—'}",f"iPTMnet count: {p.get('iptmnet_evidence_count','—')}",f"UniProt PTM B count: {p.get('uniprot_ptm_b_count','—')}"]
        self.detail.configure(state="normal"); self.detail.delete("1.0","end"); self.detail.insert("1.0","\n".join(lines)); self.detail.configure(state="disabled")

    def _export_visible_nodes(self):
        vis=self.canvas._visible_nodes()
        return [n for n in self.model["nodes"] if tuple(n["key"]) in vis]

    def _collect_export_geometry(self):
        nodes = self._export_visible_nodes()
        if not nodes:
            raise ValueError("No visible nodes to export.")
        keyset = {tuple(n["key"]) for n in nodes}
        coords = {tuple(n["key"]): self.canvas._node_world_pos(n) for n in nodes}
        edges = []
        for e in self.model["edges"]:
            a,b = tuple(e["from"]), tuple(e["to"])
            if a not in keyset or b not in keyset:
                continue
            if self.canvas.selected_path and self.canvas.selected_path not in e["path_ids"]:
                continue
            edges.append(e)
        minx = min(x for x,y in coords.values()) - 170
        maxx = max(x for x,y in coords.values()) + 170
        miny = 12
        maxy = max(y for x,y in coords.values()) + 120
        return nodes, edges, coords, keyset, minx, miny, maxx, maxy

    def _load_font(self, size, bold=False):
        if not PIL_AVAILABLE:
            return None
        candidates = []
        if os.name == "nt":
            candidates.extend([
                "C:/Windows/Fonts/arialbd.ttf" if bold else "C:/Windows/Fonts/arial.ttf",
                "C:/Windows/Fonts/segoeuib.ttf" if bold else "C:/Windows/Fonts/segoeui.ttf",
                "C:/Windows/Fonts/calibrib.ttf" if bold else "C:/Windows/Fonts/calibri.ttf",
            ])
        candidates.extend([
            "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf" if bold else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
            "/Library/Fonts/Arial Bold.ttf" if bold else "/Library/Fonts/Arial.ttf",
        ])
        for fp in candidates:
            try:
                if fp and os.path.isfile(fp):
                    return ImageFont.truetype(fp, size=size)
            except Exception:
                pass
        try:
            return ImageFont.load_default()
        except Exception:
            return None

    def _draw_arrow_head(self, draw, x1, y1, x2, y2, color, scale=1.0):
        import math
        ang = math.atan2(y2-y1, x2-x1)
        L = 11 * scale
        W = 5 * scale
        p1 = (x2, y2)
        p2 = (x2 - L*math.cos(ang) + W*math.sin(ang),
              y2 - L*math.sin(ang) - W*math.cos(ang))
        p3 = (x2 - L*math.cos(ang) - W*math.sin(ang),
              y2 - L*math.sin(ang) + W*math.cos(ang))
        draw.polygon([p1, p2, p3], fill=color)

    def _save_png_impl(self, scale=2.0):
        if not self.model:
            messagebox.showinfo("Info","Load a temporal JSON first.")
            return
        if not PIL_AVAILABLE:
            messagebox.showerror("PNG Export Error", "Pillow is required for PNG export. Please install it with: pip install pillow")
            return
        try:
            default_path = self._default_output_path(".png")
            p=filedialog.asksaveasfilename(
                title="Save PNG",
                defaultextension=".png",
                initialdir=str(Path(default_path).parent),
                initialfile=Path(default_path).name,
                filetypes=[("PNG","*.png")]
            )
            if not p:
                return

            nodes, edges, coords, keyset, minx, miny, maxx, maxy = self._collect_export_geometry()
            pad = 30
            W = int((maxx-minx)*scale + 2*pad)
            H = int((maxy-miny)*scale + 2*pad)

            img = Image.new("RGB", (max(500, W), max(400, H)), "white")
            draw = ImageDraw.Draw(img)

            # Fonts
            title_font = self._load_font(max(13, int(13*scale)), bold=True)
            gene_font = self._load_font(max(11, int(11*scale)), bold=True)
            small_font = self._load_font(max(10, int(10*scale)), bold=True)

            def tr(x, y):
                return ((x-minx)*scale + pad, (y-miny)*scale + pad)

            # Column bands and headers
            visible_cols = sorted({n["column"] for n in nodes})
            bounds=self.model.get("boundary_labels") or []
            ivs=self.model.get("intervals") or []
            for c in visible_cols:
                cx = 120 + c*self.canvas.COLUMN_GAP
                left = cx - self.canvas.COLUMN_GAP*0.46
                right = cx + self.canvas.COLUMN_GAP*0.46
                x0, y0 = tr(left, 20)
                x1, y1 = tr(right, maxy-25)
                fill = "#eff6ff" if c % 2 == 0 else "#fffaf0"
                draw.rectangle([x0,y0,x1,y1], fill=fill, outline="#e2e8f0")
                if c % 2 == 0:
                    bi=c//2; label=bounds[bi] if bi<len(bounds) else f"T{bi}"
                else:
                    ii=c//2; label=f"Bridge {ivs[ii] if ii<len(ivs) else ''}"
                tx, ty = tr(cx, 38)
                draw.text((tx, ty), label, anchor="mm", font=title_font, fill="#0f172a")

            # Edges
            for e in edges:
                a,b = tuple(e["from"]), tuple(e["to"])
                x1,y1 = coords[a]
                x2,y2 = coords[b]
                sx1, sy1 = tr(x1,y1)
                sx2, sy2 = tr(x2,y2)
                color, w = self.canvas._edge_style(e)
                width = max(1, int(round(w*scale)))
                draw.line([sx1,sy1,sx2,sy2], fill=color, width=width)
                self._draw_arrow_head(draw, sx1, sy1, sx2, sy2, color, scale)

            # Nodes
            for n in nodes:
                x,y = coords[tuple(n["key"])]
                sx, sy = tr(x,y)
                fill,border,shape = self.canvas._node_style(n)
                gene = n["gene"]
                bbox = draw.textbbox((0,0), gene, font=gene_font)
                tw = max(48*scale, (bbox[2]-bbox[0]) + 26*scale)
                th = max(28*scale, (bbox[3]-bbox[1]) + 12*scale)
                if shape == "wide_bridge":
                    bw = max(66*scale/1.0, (bbox[2]-bbox[0]) + 30*scale)
                    tip = max(10, min(bw * 0.18, 18 * scale))
                    pts = [(sx-bw/2+tip,sy-th/2),(sx+bw/2-tip,sy-th/2),(sx+bw/2,sy),(sx+bw/2-tip,sy+th/2),(sx-bw/2+tip,sy+th/2),(sx-bw/2,sy)]
                    draw.polygon(pts, fill=fill, outline=border)
                    draw.line(pts+[pts[0]], fill=border, width=max(1,int(2*scale)))
                else:
                    draw.ellipse([sx-tw/2, sy-th/2, sx+tw/2, sy+th/2], fill=fill, outline=border, width=max(1,int(2*scale)))
                txt_color = "white" if n["q_state"]=="source" else "#111827"
                draw.text((sx, sy), gene, anchor="mm", font=gene_font, fill=txt_color)


            # Export legends — large horizontal form matching the viewer.
            vs=self.model.get("validation_summary") or {}
            cell_loaded=int(vs.get("cell_records_loaded",0)); cell_matched=int(vs.get("cell_edges_matched",0))
            ph_loaded=int(vs.get("phospho_records_loaded",0)); ph_matched=int(vs.get("phospho_edges_matched",0))

            lfx,lfy=self.canvas.legend_pos["top"]
            top_x0=max(pad+4, lfx*W); top_y0=max(pad+4, lfy*H)
            top_w=min(W*0.76, 1180*scale/2.2)
            top_h=172*scale/2.2
            top_x1=min(W-pad-4, top_x0+top_w); top_y1=min(H-pad-4, top_y0+top_h)
            draw.rectangle([top_x0,top_y0,top_x1,top_y1],fill="white",outline="#94a3b8",width=max(2,int(scale)))
            draw.text((top_x0+14,top_y0+10),"NODE / Q-STATE LEGEND \n \n",font=title_font,fill="#0f172a")
            items=[("circle","#1d4ed8","#1e3a8a","Initial / boundary source"),("circle","#22c55e","#166534","Supporting anchor"),("circle","#f59e0b","#92400e","Unresolved anchor"),("circle","#8b5cf6","#5b21b6","Mixed anchor"),("circle","#ef4444","#991b1b","Opposing anchor"),("wide_bridge","#ffffff","#2563eb","Prior-knowledge bridge")]
            boxw=top_x1-top_x0-28; step=boxw/3.0
            for i,(shape,fill,border,label) in enumerate(items):
                row=i//3; col=i%3; bx=top_x0+14+col*step; by=top_y0+66+row*34*scale/2.2
                cx=bx+10*scale/2.2; cy=by+7*scale/2.2; rr=9*scale/2.2
                if shape=="wide_bridge":
                    ww=max(28*scale/2.2, rr*3.0); hh=rr*1.8; tip=max(4, rr*0.8)
                    draw.polygon([(cx-ww/2+tip,cy-hh/2),(cx+ww/2-tip,cy-hh/2),(cx+ww/2,cy),(cx+ww/2-tip,cy+hh/2),(cx-ww/2+tip,cy+hh/2),(cx-ww/2,cy)],fill=fill,outline=border)
                else: draw.ellipse([cx-rr,cy-rr,cx+rr,cy+rr],fill=fill,outline=border,width=max(1,int(scale)))
                draw.text((bx+28*scale/2.2,cy),label,anchor="lm",font=small_font,fill="#111827")

            lfx,lfy=self.canvas.legend_pos["bottom"]
            bot_x0=max(pad+4, lfx*W); bot_y0=max(pad+4, lfy*H)
            bot_w=min(W*0.84, 1330*scale/2.2)
            bot_h=220*scale/2.2
            bot_x1=min(W-pad-4, bot_x0+bot_w); bot_y1=min(H-pad-4, bot_y0+bot_h)
            draw.rectangle([bot_x0,bot_y0,bot_x1,bot_y1],fill="white",outline="#94a3b8",width=max(2,int(scale)))
            draw.text((bot_x0+14,bot_y0+10),"EDGE / VALIDATION LEGEND \n \n",font=title_font,fill="#0f172a")
            specs=[("#15803d",4,"Green: supporting or phosphosite Grade A/B"),("#d97706",4,"Orange: unresolved or phosphosite Grade C"),("#dc2626",4,"Red: opposing or cell-type UNLIKELY / NOT_EXPRESSED"),("#475569",2,"Gray-blue: general functional association")]
            boxw=bot_x1-bot_x0-28; step=boxw/2.0
            for i,(color,lw,label) in enumerate(specs):
                row=i//2; col=i%2; bx=bot_x0+14+col*step; by=bot_y0+70+row*32*scale/2.2
                draw.line([bx,by,bx+48*scale/2.2,by],fill=color,width=max(1,int(lw*scale/2.2)))
                self._draw_arrow_head(draw,bx,by,bx+48*scale/2.2,by,color,scale=max(.8,scale/2.2))
                draw.text((bx+62*scale/2.2,by),label,anchor="lm",font=small_font,fill="#111827")
            cell_txt=(f"Cell-type validation: LOADED ({cell_loaded} records; {cell_matched} matched edges)" if cell_loaded else "Cell-type validation: NOT LOADED")
            ph_txt=(f"Phosphosite validation: LOADED ({ph_loaded} records; {ph_matched} matched edges)" if ph_loaded else "Phosphosite validation: NOT LOADED")
            draw.text((bot_x0+14,bot_y1-42*scale/2.2),cell_txt,font=small_font,fill="#334155")
            draw.text((bot_x0+14,bot_y1-19*scale/2.2),ph_txt,font=small_font,fill="#334155")
            draw.text((bot_x1-14,bot_y1-19*scale/2.2),"Edge thickness ∝ repeated route occurrences",anchor="rs",font=small_font,fill="#475569")

            # Footer annotation
            info = f"Stage 5 continuous temporal viewer | {len(self.model['paths'])} paths | {len(nodes)} visible nodes | {len(edges)} visible edges"
            draw.text((pad, H-pad/2), info, anchor="ls", font=small_font, fill="#475569")
            img.save(p, format="PNG", optimize=True)
            messagebox.showinfo("Saved", f"PNG saved successfully:\n{p}")
        except Exception as e:
            messagebox.showerror("PNG Export Error", str(e))

    def save_png(self):
        self._save_png_impl(scale=2.2)

    def save_hires_png(self):
        self._save_png_impl(scale=4.0)

    def save_svg(self):
        if not self.model: messagebox.showinfo("Info","Load a temporal JSON first."); return
        try:
            default_path=self._default_output_path(".svg"); p=filedialog.asksaveasfilename(title="Save SVG",defaultextension=".svg",initialdir=str(Path(default_path).parent),initialfile=Path(default_path).name,filetypes=[("SVG","*.svg")])
            if not p:return
            nodes=self._export_visible_nodes(); ks={tuple(n["key"]) for n in nodes}; coords={tuple(n["key"]):self.canvas._node_world_pos(n) for n in nodes}
            if not nodes: raise ValueError("No visible nodes to export.")
            minx=min(x for x,y in coords.values())-140; maxx=max(x for x,y in coords.values())+140; miny=15; maxy=max(y for x,y in coords.values())+100
            parts=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{minx} {miny} {maxx-minx} {maxy-miny}">','<defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto"><path d="M0 0 L10 5 L0 10 z" fill="#475569"/></marker></defs>',f'<rect x="{minx}" y="{miny}" width="{maxx-minx}" height="{maxy-miny}" fill="white"/>']

            # export legends — large horizontal form
            vs=self.model.get("validation_summary") or {}
            cell_loaded=int(vs.get("cell_records_loaded",0)); cell_matched=int(vs.get("cell_edges_matched",0))
            ph_loaded=int(vs.get("phospho_records_loaded",0)); ph_matched=int(vs.get("phospho_edges_matched",0))
            lfx,lfy=self.canvas.legend_pos["top"]
            export_w=maxx-minx; export_h=maxy-miny
            tx0=minx + max(8, lfx*export_w); ty0=miny + max(8, lfy*export_h)
            tw=min(1180, export_w*0.76); th=128
            if tx0+tw>maxx-8: tx0=max(minx+8,maxx-8-tw)
            if ty0+th>maxy-8: ty0=max(miny+8,maxy-8-th)
            parts.append(f'<rect x="{tx0}" y="{ty0}" width="{tw}" height="{th}" fill="white" stroke="#94a3b8" stroke-width="2"/>')
            parts.append(f'<text x="{tx0+14}" y="{ty0+22}" font-family="Arial" font-size="13" font-weight="bold" fill="#0f172a">NODE / Q-STATE LEGEND</text>')
            items=[("#1d4ed8","#1e3a8a","Initial / boundary source","circle"),("#22c55e","#166534","Supporting anchor","circle"),("#f59e0b","#92400e","Unresolved anchor","circle"),("#8b5cf6","#5b21b6","Mixed anchor","circle"),("#ef4444","#991b1b","Opposing anchor","circle"),("#ffffff","#2563eb","Prior-knowledge bridge","wide_bridge")]
            step=tw/3
            for i,(fill,border,label,shape) in enumerate(items):
                row=i//3; col=i%3; bx=tx0+14+col*step; by=ty0+64+row*28
                if shape=="wide_bridge":
                    ww,hh,tip=34,18,8
                    pts=f"{bx+9-ww/2+tip},{by-hh/2} {bx+9+ww/2-tip},{by-hh/2} {bx+9+ww/2},{by} {bx+9+ww/2-tip},{by+hh/2} {bx+9-ww/2+tip},{by+hh/2} {bx+9-ww/2},{by}"
                    parts.append(f'<polygon points="{pts}" fill="{fill}" stroke="{border}" stroke-width="2"/>')
                else: parts.append(f'<circle cx="{bx+9}" cy="{by}" r="9" fill="{fill}" stroke="{border}" stroke-width="2"/>')
                parts.append(f'<text x="{bx+28}" y="{by+4}" font-family="Arial" font-size="11" font-weight="bold" fill="#111827">{xml_escape(label)}</text>')

            lfx,lfy=self.canvas.legend_pos["bottom"]
            bx0=minx + max(8, lfx*export_w); by0=miny + max(8, lfy*export_h)
            bw=min(1360, export_w*0.84); bh=144
            if bx0+bw>maxx-8: bx0=max(minx+8,maxx-8-bw)
            if by0+bh>maxy-8: by0=max(miny+8,maxy-8-bh)
            parts.append(f'<rect x="{bx0}" y="{by0}" width="{bw}" height="{bh}" fill="white" stroke="#94a3b8" stroke-width="2"/>')
            parts.append(f'<text x="{bx0+14}" y="{by0+22}" font-family="Arial" font-size="13" font-weight="bold" fill="#0f172a">EDGE / VALIDATION LEGEND</text>')
            specs=[("#15803d",4,"Green: supporting or phosphosite Grade A/B"),("#d97706",4,"Orange: unresolved or phosphosite Grade C"),("#dc2626",4,"Red: opposing or cell-type UNLIKELY / NOT_EXPRESSED"),("#475569",2,"Gray-blue: general functional association")]
            step=bw/2
            for i,(color,ew,label) in enumerate(specs):
                row=i//2; col=i%2; xx=bx0+14+col*step; yy=by0+64+row*28
                parts.append(f'<line x1="{xx}" y1="{yy}" x2="{xx+44}" y2="{yy}" stroke="{color}" stroke-width="{ew}" marker-end="url(#arrow)"/>')
                parts.append(f'<text x="{xx+58}" y="{yy+4}" font-family="Arial" font-size="11" fill="#111827">{xml_escape(label)}</text>')
            cell_txt=(f"Cell-type validation: LOADED ({cell_loaded} records; {cell_matched} matched edges)" if cell_loaded else "Cell-type validation: NOT LOADED")
            ph_txt=(f"Phosphosite validation: LOADED ({ph_loaded} records; {ph_matched} matched edges)" if ph_loaded else "Phosphosite validation: NOT LOADED")
            parts.append(f'<text x="{bx0+14}" y="{by0+126}" font-family="Arial" font-size="10" font-weight="bold" fill="#334155">{xml_escape(cell_txt)}</text>')
            parts.append(f'<text x="{bx0+14}" y="{by0+139}" font-family="Arial" font-size="10" font-weight="bold" fill="#334155">{xml_escape(ph_txt)}</text>')
            for e in self.model["edges"]:
                a,b=tuple(e["from"]),tuple(e["to"])
                if a not in ks or b not in ks:continue
                if self.canvas.selected_path and self.canvas.selected_path not in e["path_ids"]:continue
                x1,y1=coords[a]; x2,y2=coords[b]; col,w=self.canvas._edge_style(e)
                parts.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{col}" stroke-width="{w:.2f}" marker-end="url(#arrow)"/>')
            for n in nodes:
                x,y=coords[tuple(n["key"])]; fill,border,shape=self.canvas._node_style(n); gene=xml_escape(n["gene"])
                if shape=="wide_bridge":
                    ww=max(66,len(n["gene"])*8+28); hh=28; tip=max(10,min(ww*0.18,18))
                    pts=f"{x-ww/2+tip},{y-hh/2} {x+ww/2-tip},{y-hh/2} {x+ww/2},{y} {x+ww/2-tip},{y+hh/2} {x-ww/2+tip},{y+hh/2} {x-ww/2},{y}"
                    parts.append(f'<polygon points="{pts}" fill="{fill}" stroke="{border}" stroke-width="2"/>')
                else:
                    ww=max(46,len(n["gene"])*8+24); parts.append(f'<ellipse cx="{x}" cy="{y}" rx="{ww/2}" ry="14" fill="{fill}" stroke="{border}" stroke-width="2"/>')
                tc="white" if n["q_state"]=="source" else "#111827"; parts.append(f'<text x="{x}" y="{y+4}" text-anchor="middle" font-family="Arial" font-size="11" font-weight="bold" fill="{tc}">{gene}</text>')
            parts.append("</svg>"); Path(p).write_text("\n".join(parts),encoding="utf-8"); messagebox.showinfo("Saved",f"SVG saved successfully:\n{p}")
        except Exception as e: messagebox.showerror("SVG Export Error",str(e))
    def export_csv(self):
        if not self.model: messagebox.showinfo("Info","Load a temporal JSON first."); return
        try:
            default_path=self._default_output_path(".csv"); p=filedialog.asksaveasfilename(title="Save merged validation CSV",defaultextension=".csv",initialdir=str(Path(default_path).parent),initialfile=Path(default_path).name,filetypes=[("CSV","*.csv")])
            if not p:return
            fields=["interval","source","target","path_occurrences","bridge_depth","target_q_state","target_sites","cell_grade","cell_expr_a","cell_expr_b","cell_colocalization","cell_string_score","phosphosite_grade","phosphosite_confidence","matched_sites","signor_residues","signor_effects","omnipath_sources","iptmnet_count","uniprot_ptm_b_count"]
            nl={tuple(n["key"]):n for n in self.model["nodes"]}
            with open(p,"w",newline="",encoding="utf-8-sig") as f:
                w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
                for e in self.model["edges"]:
                    c=e.get("cell") or {}; ph=e.get("phospho") or {}; t=nl.get(tuple(e["to"]),{})
                    w.writerow({"interval":e["interval"],"source":e["source_gene"],"target":e["target_gene"],"path_occurrences":e.get("count",1),"bridge_depth":e.get("bridge_depth",""),"target_q_state":t.get("q_state",""),"target_sites":";".join(t.get("site_tokens") or []),"cell_grade":c.get("grade",""),"cell_expr_a":c.get("expr_a",""),"cell_expr_b":c.get("expr_b",""),"cell_colocalization":c.get("colocalization",""),"cell_string_score":c.get("string_score",""),"phosphosite_grade":ph.get("grade",""),"phosphosite_confidence":ph.get("confidence_score",""),"matched_sites":ph.get("matched_sites",""),"signor_residues":ph.get("signor_residues",""),"signor_effects":ph.get("signor_effects",""),"omnipath_sources":ph.get("omnipath_enzsub_sources") or ph.get("omnipath_interact_sources") or "","iptmnet_count":ph.get("iptmnet_evidence_count",""),"uniprot_ptm_b_count":ph.get("uniprot_ptm_b_count","")})
            messagebox.showinfo("Saved",f"CSV saved successfully:\n{p}")
        except Exception as e: messagebox.showerror("CSV Export Error",str(e))
    def open_folder(self):
        p=self.temporal_var.get().strip(); folder=os.path.dirname(p) if p else os.getcwd()
        try:
            if os.name=="nt": os.startfile(folder)
            elif sys.platform=="darwin":
                import subprocess; subprocess.Popen(["open",folder])
            else:
                import subprocess; subprocess.Popen(["xdg-open",folder])
        except Exception as e: messagebox.showerror("Open Folder Error",str(e))

if __name__=="__main__":
    App().mainloop()
