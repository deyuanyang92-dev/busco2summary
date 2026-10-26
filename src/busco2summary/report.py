"""Self-contained interactive HTML report (no network, no CDN)."""

from __future__ import annotations

import html
import json
import time
from pathlib import Path
from typing import List

from . import __version__
from .analysis import Dataset

_CSS = r"""
:root{--bg:#fcfcfb;--panel:#ffffff;--ink:#0b0b0b;--ink2:#52514e;--muted:#8a8984;--line:#e4e3df;
--accent:#2a78d6;--S:#4A9FD8;--D:#1F5F99;--F:#D08F00;--M:#C93C3C;--warn:#b25f00;--bad:#b3261e;--ok:#1b7f3b;
--setS:#2a78d6;--setSD:#6250d6;--setSF:#eb6834;--setSDF:#1baf7a;color-scheme:light}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#1a1a19;--panel:#222221;--ink:#f0efec;
--ink2:#c3c2b7;--muted:#8f8e86;--line:#3a3a37;--accent:#3987e5;--S:#3F95D3;--D:#2D62A6;--F:#BA7E00;--M:#C8405A;
--warn:#e0a040;--bad:#ef6a5f;--ok:#4cbf6b;--setS:#3987e5;--setSD:#9085e9;--setSF:#d95926;--setSDF:#199e70;color-scheme:dark}}
:root[data-theme="dark"]{--bg:#1a1a19;--panel:#222221;--ink:#f0efec;--ink2:#c3c2b7;--muted:#8f8e86;--line:#3a3a37;
--accent:#3987e5;--S:#3F95D3;--D:#2D62A6;--F:#BA7E00;--M:#C8405A;--warn:#e0a040;--bad:#ef6a5f;--ok:#4cbf6b;
--setS:#3987e5;--setSD:#9085e9;--setSF:#d95926;--setSDF:#199e70;color-scheme:dark}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);
font:14px/1.45 system-ui,-apple-system,"Segoe UI",Roboto,"Noto Sans SC",sans-serif}
header{padding:20px 24px 8px;max-width:1280px;margin:auto}h1{font-size:20px;margin:0 0 4px}
.sub{color:var(--ink2);font-size:13px}nav{display:flex;gap:4px;flex-wrap:wrap;padding:8px 24px;max-width:1280px;margin:auto;
border-bottom:1px solid var(--line)}nav button{background:none;border:0;padding:8px 12px;border-radius:6px;color:var(--ink2);
font:inherit;cursor:pointer}nav button.on{background:var(--panel);color:var(--ink);box-shadow:0 0 0 1px var(--line);font-weight:600}
main{max-width:1280px;margin:auto;padding:16px 24px 48px}section{display:none}section.on{display:block}
.cards{display:grid;grid-template-columns:repeat(auto-fill,minmax(150px,1fr));gap:10px;margin-bottom:16px}
.card{background:var(--panel);border:1px solid var(--line);border-radius:8px;padding:10px 12px}
.card .v{font-size:22px;font-weight:650;font-variant-numeric:tabular-nums}.card .k{color:var(--ink2);font-size:12px}
.panel{background:var(--panel);border:1px solid var(--line);border-radius:8px;padding:14px;margin-bottom:16px;overflow:auto}
h2{font-size:15px;margin:0 0 10px}table{border-collapse:collapse;width:100%;font-size:12.5px;font-variant-numeric:tabular-nums}
th,td{padding:4px 8px;border-bottom:1px solid var(--line);text-align:right;white-space:nowrap}th{position:sticky;top:0;
background:var(--panel);cursor:pointer;color:var(--ink2);font-weight:600}td:first-child,th:first-child{text-align:left}
tr:hover td{background:color-mix(in srgb,var(--accent) 7%,transparent)}.tbl{max-height:620px;overflow:auto}
input[type=search],select{background:var(--panel);color:var(--ink);border:1px solid var(--line);border-radius:6px;padding:6px 8px;font:inherit}
.row{display:flex;gap:16px;flex-wrap:wrap;align-items:center;margin-bottom:10px}.row label{color:var(--ink2)}
input[type=range]{width:220px;accent-color:var(--accent)}.btn{background:var(--accent);color:#fff;border:0;border-radius:6px;
padding:7px 12px;font:inherit;cursor:pointer}.btn.ghost{background:none;color:var(--accent);box-shadow:inset 0 0 0 1px var(--accent)}
.legend{display:flex;gap:14px;flex-wrap:wrap;font-size:12px;color:var(--ink2);margin:6px 0}.sw{display:inline-block;width:10px;height:10px;
border-radius:2px;margin-right:5px;vertical-align:-1px}.flag{display:inline-block;padding:0 6px;border-radius:9px;font-size:11px;
margin-left:3px;border:1px solid var(--line);color:var(--ink2)}.flag.bad{color:var(--bad);border-color:var(--bad)}
.flag.warn{color:var(--warn);border-color:var(--warn)}#tip{position:fixed;pointer-events:none;background:var(--panel);
border:1px solid var(--line);border-radius:6px;padding:6px 8px;font-size:12px;box-shadow:0 4px 14px rgba(0,0,0,.15);display:none;z-index:9}
pre{white-space:pre-wrap;color:var(--ink2);font-size:12.5px;margin:0}.note{color:var(--ink2);font-size:12.5px}
svg text{fill:var(--ink2);font-size:11px}.big{font-size:28px;font-weight:700}
@media (max-width:640px){header,nav,main{padding-left:16px;padding-right:16px}input[type=range]{width:150px}}
"""

_JS = r"""
const D=JSON.parse(document.getElementById('data').textContent);
const $=s=>document.querySelector(s), tip=$('#tip');
const N=D.samples.length, G=D.genes.length, ST=['S','D','F','M'];
const fmtP=x=>(x*100).toFixed(1)+'%';
function show(id){document.querySelectorAll('nav button').forEach(b=>b.classList.toggle('on',b.dataset.t===id));
 document.querySelectorAll('section').forEach(s=>s.classList.toggle('on',s.id===id)); if(id==='mx')drawMatrix(); if(id==='occ')updOcc();}
document.querySelectorAll('nav button').forEach(b=>b.onclick=()=>show(b.dataset.t));
function tipAt(e,h){tip.innerHTML=h;tip.style.display='block';const x=Math.min(e.clientX+14,innerWidth-260);tip.style.left=x+'px';tip.style.top=(e.clientY+14)+'px'}
function tipOff(){tip.style.display='none'}
function dl(name,text){const a=document.createElement('a');a.href=URL.createObjectURL(new Blob([text],{type:'text/plain'}));a.download=name;a.click()}
function sortable(tbl,rows,cols,render){let key=null,asc=true;const th=tbl.querySelector('thead tr');
 th.innerHTML=cols.map(c=>`<th data-k="${c[0]}">${c[1]}</th>`).join('');
 th.querySelectorAll('th').forEach(h=>h.onclick=()=>{const k=h.dataset.k;asc=key===k?!asc:false;key=k;
  rows.sort((a,b)=>{const x=a[k],y=b[k];if(x===y)return 0;if(x==null)return 1;if(y==null)return-1;return (x<y?-1:1)*(asc?1:-1)});render()});}
// ---------- per gene counts ----------
const cnt=D.rows.map(r=>{const c={S:0,D:0,F:0,M:0};for(const ch of r)c[ch]++;return c});
// ---------- overview ----------
(function(){const st={};D.samples.forEach(s=>st[s.state]=(st[s.state]||0)+1);
 const mean=k=>D.samples.reduce((a,s)=>a+s['pct_'+k],0)/Math.max(N,1);
 const cards=[[N,'samples'],[G,'BUSCOs ('+D.meta.lineage+')'],[fmtP(mean('C')),'mean complete (C)'],[fmtP(mean('S')),'mean single-copy (S)'],
  [fmtP(mean('D')),'mean duplicated (D)'],[fmtP(mean('F')),'mean fragmented (F)'],[fmtP(mean('M')),'mean missing (M)'],
  [cnt.filter(c=>c.S===N).length,'BUSCOs S in all samples'],[cnt.filter(c=>c.S+c.F===N).length,'BUSCOs S or F in all']];
 $('#cards').innerHTML=cards.map(c=>`<div class="card"><div class="v">${c[0]}</div><div class="k">${c[1]}</div></div>`).join('');
 $('#report').textContent=D.meta.report.join('\n');
 const fl=D.samples.filter(s=>s.qc_flags);$('#qcn').textContent=fl.length;
 $('#qc').innerHTML='<table><thead><tr><th>sample</th><th>flags</th><th>C</th><th>F</th><th>M</th></tr></thead><tbody>'+
  fl.map(s=>`<tr><td>${s.sample}</td><td>${flags(s.qc_flags)}</td><td>${fmtP(s.pct_C)}</td><td>${fmtP(s.pct_F)}</td><td>${fmtP(s.pct_M)}</td></tr>`).join('')+'</tbody></table>';})();
function flags(f){return (f||'').split(',').filter(Boolean).map(x=>`<span class="flag ${/EMPTY|FAIL|LOW_C|MISMATCH/.test(x)?'bad':'warn'}">${x}</span>`).join('')}
// ---------- samples ----------
(function(){const rows=D.samples.slice();const tbl=$('#stbl');
 const cols=[['sample','sample'],['bar','S / D / F / M'],['pct_C','C'],['pct_S','S'],['pct_D','D'],['pct_F','F'],['pct_M','M'],['pct_SF','S+F'],
  ['mean_F_ratio','F len ratio'],['pct_E','E %'],['mode','mode'],['busco_version','BUSCO'],['qc_flags','QC']];
 function bar(s){let x=0;return '<svg width="160" height="12" role="img" aria-label="S '+fmtP(s.pct_S)+'">'+ST.map(k=>{const w=s['pct_'+k]*160;const r=`<rect x="${x}" y="1" width="${Math.max(w-1,0)}" height="10" fill="var(--${k})"/>`;x+=w;return r}).join('')+'</svg>'}
 function render(){const q=$('#sq').value.toLowerCase();tbl.querySelector('tbody').innerHTML=rows.filter(s=>!q||s.sample.toLowerCase().includes(q)||(s.qc_flags||'').toLowerCase().includes(q)).map(s=>
  `<tr><td>${s.sample}</td><td>${bar(s)}</td><td>${fmtP(s.pct_C)}</td><td>${fmtP(s.pct_S)}</td><td>${fmtP(s.pct_D)}</td><td>${fmtP(s.pct_F)}</td><td>${fmtP(s.pct_M)}</td><td>${fmtP(s.pct_SF)}</td><td>${s.mean_F_ratio==null?'':s.mean_F_ratio.toFixed(2)}</td><td>${s.pct_E==null?'':s.pct_E}</td><td>${s.mode||''}</td><td>${s.busco_version||''}</td><td>${flags(s.qc_flags)}</td></tr>`).join('')}
 sortable(tbl,rows,cols,render);$('#sq').oninput=render;render();})();
// ---------- genes ----------
(function(){const rows=D.genes.map((g,i)=>({id:g,S:cnt[i].S,D:cnt[i].D,F:cnt[i].F,M:cnt[i].M,occS:cnt[i].S/N,occSF:(cnt[i].S+cnt[i].F)/N,
 occSDF:(N-cnt[i].M)/N,gain:cnt[i].F/N,fr:D.gfr[i],fl:D.gflags[i]}));const tbl=$('#gtbl');
 const cols=[['id','BUSCO'],['S','n S'],['D','n D'],['F','n F'],['M','n M'],['occS','occ S'],['occSF','occ S+F'],['occSDF','occ any'],['gain','F gain'],['fr','median F ratio'],['fl','flags']];
 function render(){const q=$('#gq').value.toLowerCase();const f=rows.filter(r=>!q||r.id.toLowerCase().includes(q)||r.fl.toLowerCase().includes(q));
  $('#gcount').textContent=f.length+' BUSCOs'+(f.length>3000?' (first 3000 shown)':'');
  tbl.querySelector('tbody').innerHTML=f.slice(0,3000).map(r=>`<tr><td>${r.id}</td><td>${r.S}</td><td>${r.D}</td><td>${r.F}</td><td>${r.M}</td><td>${fmtP(r.occS)}</td><td>${fmtP(r.occSF)}</td><td>${fmtP(r.occSDF)}</td><td>${fmtP(r.gain)}</td><td>${r.fr==null?'':r.fr.toFixed(2)}</td><td>${flags(r.fl)}</td></tr>`).join('')}
 sortable(tbl,rows,cols,render);$('#gq').oninput=render;render();})();
// ---------- occupancy explorer ----------
const SETS={S:'S',SD:'SD',SF:'SF',SDF:'SDF'};
function present(r,j,set){return set.includes(r[j])}
function occCompute(set,minComp,minOcc){const comp=D.samples.map((s,j)=>{let c=0;for(const r of D.rows)if(present(r,j,set))c++;return c/Math.max(G,1)});
 const keep=[];comp.forEach((c,j)=>{if(c>=minComp)keep.push(j)});const n=keep.length;const need=Math.ceil(minOcc*n-1e-9);
 const genes=[];let cells=0;const per=D.rows.map(r=>{let c=0;for(const j of keep)if(present(r,j,set))c++;return c});
 per.forEach((c,i)=>{if(c>=need&&c>0){genes.push(D.genes[i]);cells+=c}});
 return {keep,n,genes,fill:genes.length&&n?cells/(genes.length*n):0,per,comp}}
let OC=null;
function updOcc(){const set=$('#oset').value,mc=+$('#ocomp').value/100,mo=+$('#oocc').value/100;
 $('#ocompv').textContent=Math.round(mc*100)+'%';$('#ooccv').textContent=Math.round(mo*100)+'%';
 OC=occCompute(set,mc,mo);$('#ores').innerHTML=`<div class="cards"><div class="card"><div class="v">${OC.n}</div><div class="k">samples kept (of ${N})</div></div>
 <div class="card"><div class="v">${OC.genes.length}</div><div class="k">BUSCOs in ≥${Math.ceil(mo*OC.n-1e-9)} samples</div></div>
 <div class="card"><div class="v">${fmtP(OC.fill)}</div><div class="k">matrix fill (present cells)</div></div></div>`;
 // threshold table for all sets
 const ts=[.5,.6,.7,.75,.8,.85,.9,.95,1];let h='<table><thead><tr><th>occupancy ≥</th>'+Object.keys(SETS).map(k=>`<th>${k}</th>`).join('')+'<th>F gain (SF−S)</th></tr></thead><tbody>';
 const per={};for(const k in SETS){per[k]=D.rows.map(r=>{let c=0;for(const j of OC.keep)if(present(r,j,SETS[k]))c++;return c})}
 for(const t of ts){const need=Math.ceil(t*OC.n-1e-9);const v={};for(const k in SETS)v[k]=per[k].filter(c=>c>=need&&c>0).length;
  h+=`<tr><td>${Math.round(t*100)}% (≥${need})</td>`+Object.keys(SETS).map(k=>`<td>${v[k]}</td>`).join('')+`<td>+${v.SF-v.S}</td></tr>`}
 $('#otbl').innerHTML=h+'</tbody></table>';drawCurve(per);}
function drawCurve(per){const W=Math.min($('#ocurve').clientWidth||800,1000),H=300,L=50,R=40,T=10,B=36;const n=OC.n||1;
 const x=v=>L+(W-L-R)*v,y=v=>T+(H-T-B)*(1-v/Math.max(G,1));let s=`<svg width="${W}" height="${H}">`;
 for(let g=0;g<=4;g++){const v=G*g/4;s+=`<line x1="${L}" x2="${W-R}" y1="${y(v)}" y2="${y(v)}" stroke="var(--line)"/><text x="${L-6}" y="${y(v)+4}" text-anchor="end">${Math.round(v)}</text>`}
 for(let t=0;t<=10;t++){s+=`<text x="${x(t/10)}" y="${H-B+16}" text-anchor="middle">${t*10}%</text>`}
 s+=`<text x="${(L+W-R)/2}" y="${H-4}" text-anchor="middle">minimum % of kept samples containing the BUSCO</text>`;
 for(const k in SETS){const hist=new Array(n+1).fill(0);per[k].forEach(c=>hist[c]++);let acc=0;const pts=[];
  for(let c=n;c>=1;c--){acc+=hist[c];pts.push([c/n,acc])}pts.reverse();
  s+=`<polyline fill="none" stroke="var(--set${k})" stroke-width="2" points="${pts.map(p=>x(p[0])+','+y(p[1])).join(' ')}"/>`;
  const half=pts.find(p=>p[1]<=G/2);if(half)s+=`<text x="${x(half[0])+4}" y="${y(G/2)-4}">${k}</text>`}
 const mo=+$('#oocc').value/100;s+=`<line x1="${x(mo)}" x2="${x(mo)}" y1="${T}" y2="${H-B}" stroke="var(--ink2)" stroke-dasharray="3 3"/>`;
 s+=`<rect x="${L}" y="${T}" width="${W-L-R}" height="${H-T-B}" fill="transparent" id="ohit"/></svg>`;$('#ocurve').innerHTML=s;
 $('#ohit').onmousemove=e=>{const r=e.target.getBoundingClientRect();const v=Math.max(0,Math.min(1,(e.clientX-r.left)/r.width));const need=Math.ceil(v*n-1e-9);
  tipAt(e,`<b>≥${Math.round(v*100)}%</b> (≥${need} samples)<br>`+Object.keys(SETS).map(k=>`${k}: ${per[k].filter(c=>c>=need&&c>0).length}`).join('<br>'))};
 $('#ohit').onmouseleave=tipOff;}
['#oset','#ocomp','#oocc'].forEach(s=>$(s).oninput=updOcc);
$('#odl').onclick=()=>dl(`busco_${$('#oset').value}_occ${$('#oocc').value}_comp${$('#ocomp').value}.txt`,OC.genes.join('\n')+'\n');
$('#odls').onclick=()=>dl(`samples_kept_comp${$('#ocomp').value}.txt`,OC.keep.map(j=>D.samples[j].sample).join('\n')+'\n');
// ---------- matrix ----------
function drawMatrix(){const cv=$('#mcv');
 const so=D.samples.map((s,j)=>j).sort((a,b)=>D.samples[b].pct_SDF-D.samples[a].pct_SDF);
 const go=D.genes.map((g,i)=>i).sort((a,b)=>(cnt[a].M-cnt[b].M)||(cnt[a].F-cnt[b].F));
 const cw=Math.max(2,Math.min(14,Math.floor(1100/Math.max(N,1)))),ch=Math.max(1,Math.min(6,Math.floor(900/Math.max(G,1))));
 cv.width=N*cw;cv.height=G*ch;const ctx=cv.getContext('2d');const cs=getComputedStyle(document.documentElement);
 const col={S:cs.getPropertyValue('--S'),D:cs.getPropertyValue('--D'),F:cs.getPropertyValue('--F'),M:cs.getPropertyValue('--M')};
 go.forEach((gi,yy)=>{const r=D.rows[gi];so.forEach((j,xx)=>{ctx.fillStyle=col[r[j]];ctx.fillRect(xx*cw,yy*ch,cw-(cw>3?1:0),ch)})});
 cv.onmousemove=e=>{const b=cv.getBoundingClientRect();const xx=Math.floor((e.clientX-b.left)/b.width*N),yy=Math.floor((e.clientY-b.top)/b.height*G);
  if(xx<0||yy<0||xx>=N||yy>=G)return tipOff();const j=so[xx],gi=go[yy];tipAt(e,`<b>${D.genes[gi]}</b><br>${D.samples[j].sample}<br>status: <b>${D.rows[gi][j]}</b>`)};
 cv.onmouseleave=tipOff;}
$('#theme').onclick=()=>{const r=document.documentElement;const dark=r.dataset.theme?r.dataset.theme==='dark':matchMedia('(prefers-color-scheme: dark)').matches;
 r.dataset.theme=dark?'light':'dark';if($('#mx').classList.contains('on'))drawMatrix();};
show('ov');
"""


def write_report(ds: Dataset, sample_rows: List[dict], gene_rows: List[dict],
                 report_lines: List[str], out: Path) -> Path:
    keep = ("sample", "state", "pct_S", "pct_D", "pct_C", "pct_F", "pct_M", "pct_SF",
            "pct_SDF", "mean_F_ratio", "qc_flags")
    samples = []
    for s, r in zip(ds.samples, sample_rows):
        d = {k: r.get(k) for k in keep}
        d.update(pct_E=s.info.get("pct_E"), mode=s.info.get("mode"),
                 busco_version=s.info.get("busco_version") or s.tool)
        samples.append(d)
    data = {
        "meta": {"version": __version__, "lineage": ds.lineage, "ref": ds.ref_len_source,
                 "created": time.strftime("%Y-%m-%d %H:%M"), "report": report_lines},
        "samples": samples,
        "genes": ds.genes,
        "rows": [ds.rows[g] for g in ds.genes],
        "gfr": [None if r["median_F_ratio"] is None else round(r["median_F_ratio"], 3)
                for r in gene_rows],
        "gflags": [r["flags"] for r in gene_rows],
    }
    blob = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    title = f"BUSCO summary - {html.escape(ds.lineage)}"
    page = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{title}</title><style>{_CSS}</style></head><body>
<header><div class="row" style="justify-content:space-between;margin:0"><h1>{title}</h1>
<button class="btn ghost" id="theme" title="toggle light/dark">◐ theme</button></div>
<div class="sub">busco2summary v{__version__} · {len(ds.samples)} samples · {len(ds.genes)} BUSCOs ·
fragment reference length: {html.escape(ds.ref_len_source)} · {time.strftime('%Y-%m-%d %H:%M')}</div></header>
<nav><button data-t="ov">Overview</button><button data-t="sm">Samples</button><button data-t="gn">Genes (S/D/F/M)</button>
<button data-t="occ">Occupancy &amp; gene lists</button><button data-t="mx">Matrix</button></nav>
<main>
<section id="ov"><div class="cards" id="cards"></div>
<div class="panel"><h2>Dataset</h2><pre id="report"></pre></div>
<div class="panel"><h2>Samples with QC flags (<span id="qcn"></span>)</h2><div class="tbl" id="qc"></div>
<p class="note">LOW_C complete &lt; threshold · HIGH_D duplicated high (haplotigs, contamination, isoforms) · HIGH_F many fragments
(short contigs, low coverage) · HIGH_M · HIGH_E internal stop codons (BUSCO ≥ 5.8) · EMPTY no BUSCO found · MISMATCH full_table ≠ short_summary.</p></div></section>
<section id="sm"><div class="panel"><div class="row"><input type="search" id="sq" placeholder="filter samples or flags…">
<div class="legend"><span><i class="sw" style="background:var(--S)"></i>S single-copy</span><span><i class="sw" style="background:var(--D)"></i>D duplicated</span>
<span><i class="sw" style="background:var(--F)"></i>F fragmented</span><span><i class="sw" style="background:var(--M)"></i>M missing</span></div></div>
<div class="tbl"><table id="stbl"><thead><tr></tr></thead><tbody></tbody></table></div>
<p class="note">Click a column header to sort. "F len ratio" = mean length of fragmented BUSCOs relative to the reference length.</p></div></section>
<section id="gn"><div class="panel"><div class="row"><input type="search" id="gq" placeholder="BUSCO id or flag (F_PRONE, D_PRONE…)"><span class="note" id="gcount"></span></div>
<div class="tbl"><table id="gtbl"><thead><tr></tr></thead><tbody></tbody></table></div>
<p class="note">F gain = fraction of samples where the BUSCO is only fragmented: occupancy you recover by accepting fragments.
F_PRONE ≥30 % samples fragmented · D_PRONE ≥20 % duplicated (paralogy risk) · F_RESCUED_50 passes 50 % occupancy only with F.</p></div></section>
<section id="occ"><div class="panel"><h2>Occupancy explorer</h2>
<div class="row"><label>status set <select id="oset"><option value="S">S - single-copy</option><option value="SD">SD - complete</option>
<option value="SF" selected>SF - single + fragmented</option><option value="SDF">SDF - any hit</option></select></label>
<label>min sample completeness <input type="range" id="ocomp" min="0" max="100" value="0"> <b id="ocompv"></b></label>
<label>min gene occupancy <input type="range" id="oocc" min="0" max="100" value="50"> <b id="ooccv"></b></label></div>
<div id="ores"></div><div class="row"><button class="btn" id="odl">Download BUSCO IDs</button><button class="btn ghost" id="odls">Download kept samples</button></div>
<div class="legend"><span><i class="sw" style="background:var(--setS)"></i>S single-copy</span><span><i class="sw" style="background:var(--setSD)"></i>SD complete</span>
<span><i class="sw" style="background:var(--setSF)"></i>SF single + fragmented</span><span><i class="sw" style="background:var(--setSDF)"></i>SDF any hit</span></div>
<div id="ocurve" style="width:100%"></div><div class="tbl" id="otbl"></div>
<p class="note">A BUSCO passes when present in at least ceil(occupancy × kept samples) samples. Samples are first filtered by their
completeness under the chosen status set. The same lists are written to <code>gene_lists/</code> by the command line.</p></div></section>
<section id="mx"><div class="panel"><div class="row"><span class="note">Columns: samples (most complete left) · rows: BUSCOs (fewest missing on top). Hover a cell.</span>
<div class="legend"><span><i class="sw" style="background:var(--S)"></i>S</span><span><i class="sw" style="background:var(--D)"></i>D</span>
<span><i class="sw" style="background:var(--F)"></i>F</span><span><i class="sw" style="background:var(--M)"></i>M</span></div></div>
<div style="overflow:auto;max-height:80vh"><canvas id="mcv" style="width:100%;image-rendering:pixelated"></canvas></div></div></section>
</main><div id="tip"></div>
<script type="application/json" id="data">{blob}</script><script>{_JS}</script></body></html>"""
    out.write_text(page, encoding="utf-8")
    return out
