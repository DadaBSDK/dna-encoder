"""Build a standalone LaTeX report from saved example and benchmark evidence."""
from pathlib import Path
import csv,json,hashlib

ROOT=Path(__file__).resolve().parent.parent
EX=ROOT/'results/report_examples_20261008'
VALID=ROOT/'results/phase4/validation_optimized_20261008'
x=json.loads((EX/'examples.json').read_text())
v=list(csv.DictReader((VALID/'summary.csv').open()))
scale=list(csv.DictReader((ROOT/'results/phase3/cluster_scale_100KB/cluster_scale.csv').open()))


def esc(s):
    return str(s).replace('\\',r'\textbackslash{}').replace('_',r'\_').replace('%',r'\%').replace('&',r'\&').replace('#',r'\#')

arms=['whiten-RS','goldman','steering-A-P8','steering-C-P8','steering-B-P6','fountain']
names=['Whiten + RS','Goldman','Steering A, P=8','Steering C, P=8','Steering B, P=6','Fountain']
lookup={(r['arm'],r['payload_seed'],r['strength']):r for r in v}
vt=[]
for a,n in zip(arms,names):
    vt.append(n+' & '+' & '.join(lookup[a,str(p),s]['recovered']+'/20' for p in (0,1) for s in ('erasure','repair'))+r' \\')
st=[]
for a,n in zip(['naive2bit-whiten',*arms[1:]],names):
    b=next(r for r in scale if r['arm']==a and r['channel']=='badread')
    ids=next(r for r in scale if r['arm']==a and r['channel']=='ids9')
    st.append(f"{n} & {int(b['n_oligos']):,} & {float(b['min_ed_frac']):.4f} & {100*float(ids['split_oligo_frac']):.2f}"+r' \\')
examples={r['example']:r for r in x['examples'] if r['strength']=='repair'}
rows=[]
labels=['Text, noiseless','Binary, noiseless','Two oligos removed','3\\% IDS, depth 12','9\\% IDS, depth 2']
for i,label in enumerate(labels,1):
    r=examples[f'E{i}']
    result=r'\textcolor{teal}{Pass}' if r['recovered'] else r'\textcolor{red}{Fail}'
    rows.append(f"E{i} & {label} & {r['input_bytes']:,} & {r['reads']:,} & {result} & {result}"+r' \\')
errors=[]
for k in ('substitutions','insertions','deletions'):
    errors.append(' '.join(f"({eid},{100*examples[eid][k]/(examples[eid]['reads']*200):.5f})" for eid in ('E4','E5')))
dna=x['pools']['text']['first_data_oligo']
hashes='\n'.join(name+':\n'+x['inputs'][name]['sha256'][:32]+'\n'+x['inputs'][name]['sha256'][32:] for name in ('text','binary'))

# ---- Phase 2 (historical, format v2 before the inner CRC-32 change)
import ast,re
import xml.etree.ElementTree as ET
SWEEP=ROOT/'results/phase2/sweep_v2/summary.csv'
MFE=ROOT/'results/phase2/composition/mfe_summary.csv'
sweep={r['arm']:r for r in csv.DictReader(SWEEP.open())}
p2=[]
for a,n in [('naive2bit','Quaternary, raw'),('naive2bit-whiten','Quaternary, whitened'),('goldman','Goldman rotating ternary'),
            ('steering-A-P8','Steering A, P=8'),('steering-C-P8','Steering C, P=8'),('steering-B-P6','Steering B, P=6')]:
    r=sweep[a]
    p2.append(f"{n} & {float(r['code_density']):.3f} & {float(r['effective_density']):.3f} & {r['max_hp']} & {100*float(r['frac_violating']):.1f}\\%"+r' \\')
mfe={(r['source'],float(r['temp'])):r for r in csv.DictReader(MFE.open())}
def mfe_cell(src,t):
    r=mfe[src,t]
    return f"{float(r['mfe_median']):.1f} [{float(r['mfe_median_ci_lo']):.1f}, {float(r['mfe_median_ci_hi']):.1f}]"
mfe_rows=[f"{label} & {mfe_cell(src,37.0)} & {mfe_cell(src,60.0)} & {float(mfe[src,37.0]['paired_frac_mean']):.2f}"+r' \\'
          for src,label in [('iid','Whitened quaternary (i.i.d.)'),('rot','Rotating ternary (run length 1)'),('shiftban','Control: different banned successor')]]

def md_row(path,label):
    for line in path.read_text().splitlines():
        if line.startswith('|') and label in line:
            return [c.strip().replace('**','') for c in line.strip().strip('|').split('|')]
    raise KeyError(label)
STRUCT=ROOT/'docs/notes/structure.md'
def rho(label):
    cell=md_row(STRUCT,label)[2]
    return re.search(r'[+\u2212-]\d\.\d+',cell).group().replace('\u2212','-')
proxy=[('Incremental k-mer proxy (tier 1)',*rho('Tier 1, incremental k-mer').split()),
       ('Local partition function, window 80 nt (tier 2)',rho('Tier 2, local pf, W = 80')),
       ('Global MFE at 60 \\textdegree C',rho('Global MFE min(fwd, rc) at 60'))]
proxy_rows=[f"{n} & ${v.replace('-','-').replace('+','+')}$"+r' \\' for n,v in proxy]
CAL=ROOT/'docs/notes/channel_calibration.md'
def cal(label):
    cells=md_row(CAL,label)
    return ' & '.join(esc(c).replace('—','--') for c in cells[1:])
cal_rows=[f"Per-base error & {cal('Per-base error')}"+r' \\',
          f"Sub.\\,/\\,del.\\,/\\,ins. & {cal('Substitution / deletion / insertion')}"+r' \\']

# ---- Testing campaign (rerun for this report)
TEST=ROOT/'results/testing_20261009'
cases=list(ET.parse(TEST/'pytest_junit.xml').getroot().iter('testcase'))
n_xfail=sum(1 for c in cases for k in c.findall('skipped') if 'xfail' in (k.get('type','')+k.get('message','')))
n_skip=sum(1 for c in cases if c.find('skipped') is not None)-n_xfail
n_fail=sum(1 for c in cases if c.find('failure') is not None or c.find('error') is not None)
n_pass=len(cases)-n_xfail-n_skip-n_fail
coverage=json.loads((TEST/'coverage.json').read_text())['totals']['percent_covered']
tap=(TEST/'node_tests.tap').read_text()
node_pass=int(re.search(r'^# pass (\d+)',tap,re.M).group(1)); node_fail=int(re.search(r'^# fail (\d+)',tap,re.M).group(1))
ADV=ROOT/'tests/test_adversarial.py'
bugs=[]
for node in ast.walk(ast.parse(ADV.read_text())):
    if isinstance(node,ast.FunctionDef):
        for d in node.decorator_list:
            if isinstance(d,ast.Call) and getattr(d.func,'attr','')=='xfail':
                reason=next(k.value.value for k in d.keywords if k.arg=='reason')
                bugs.append((node.name,reason.removeprefix('BUG: ')))
guards=sum(1 for node in ast.walk(ast.parse(ADV.read_text())) if isinstance(node,ast.FunctionDef) and node.name.startswith('test_')
           and not any(isinstance(d,ast.Call) and getattr(d.func,'attr','')=='xfail' for d in node.decorator_list))
def tex_text(t):
    return re.sub(r'`([^`]*)`',lambda m:r'\texttt{'+m.group(1)+'}',esc(t).replace('~',r'\textasciitilde{}'))
bug_rows=[f"B{i} & {tex_text(r)}"+r' \\' for i,(_,r) in enumerate(bugs,1)]

# ---- Web demo: error-tolerant browser codec benchmark
SITE=ROOT/'results/site_resilience_20261009'
site=list(csv.DictReader((SITE/'summary.csv').open()))
site_prov=json.loads((SITE/'provenance.json').read_text())
site_trials=sum(int(r['trials']) for r in site)
S={(int(r['parity_permille']),float(r['error_per_base']),int(r['reads_per_strand'])):r for r in site}
errs=sorted({float(r['error_per_base']) for r in site}); readsv=sorted({int(r['reads_per_strand']) for r in site})
parities=sorted({int(r['parity_permille']) for r in site})
site_grid=[]
for e in errs:
    for rd in readsv:
        cells=' & '.join(f"{S[p,e,rd]['exact']}/{S[p,e,rd]['trials']}" for p in parities)
        site_grid.append(f"{100*e:g}\\% & {rd} & {cells}"+r' \\')
site_partial=[]
for e in errs:
    for rd in readsv:
        r=S[250,e,rd]
        site_partial.append(f"{100*e:g}\\% & {rd} & {r['exact']}/{r['trials']} & {100*float(r['mean_guaranteed_exact']):.1f}\\% & {100*float(r['mean_actually_correct']):.1f}\\% & {float(r['median_decode_ms']):.0f}"+r' \\')
site_plot='\n'.join(r'\addplot+[mark=*,thick] coordinates {'+' '.join(f"({rd},{S[250,e,rd]['exact']})" for rd in readsv)+'};'
                    +r'\addlegendentry{'+f'{100*e:g}'+r'\% per base}' for e in errs)
nt_per_byte={p:float(S[p,errs[0],readsv[0]]['nt_per_byte']) for p in parities}
real=S[250,0.01,3]

tex=r'''\documentclass[11pt,a4paper]{article}
\usepackage[T1]{fontenc}
\usepackage{lmodern,microtype}
\usepackage[margin=22mm,headheight=15pt]{geometry}
\usepackage{amsmath,amssymb,booktabs,tabularx,array,xcolor}
\usepackage{tikz,pgfplots}
\pgfplotsset{compat=1.18}
\usetikzlibrary{arrows.meta,positioning}
\usepackage{fancyhdr,hyperref,listings}
\definecolor{navy}{HTML}{17324D}
\definecolor{teal}{HTML}{087E8B}
\definecolor{light}{HTML}{EFF5F8}
\definecolor{red}{HTML}{B33C35}
\definecolor{muted}{HTML}{526577}
\hypersetup{colorlinks=true,linkcolor=teal,urlcolor=teal,pdftitle={DNA Data Storage Encoder: Project Report},pdfauthor={}}
\pagestyle{fancy}\fancyhf{}
\fancyhead[L]{\small\textcolor{navy}{DNA DATA STORAGE ENCODER}}
\fancyhead[R]{\small\textcolor{muted}{Project report | 09 Oct 2026}}
\fancyfoot[L]{\footnotesize Implementation v0.1.0 / format v3}
\fancyfoot[R]{\footnotesize\thepage}
\renewcommand{\headrulewidth}{0.4pt}
\setlength{\parindent}{0pt}\setlength{\parskip}{6pt}
\setlength{\tabcolsep}{6pt}\renewcommand{\arraystretch}{1.18}
\newcommand{\note}[1]{\par\noindent\colorbox{light}{\parbox{\dimexpr\linewidth-2\fboxsep\relax}{\small #1}}\par}
\newcommand{\pageheading}[2]{\clearpage\section{#1}{\color{muted}\small #2}\par\vspace{5pt}}
\lstset{basicstyle=\ttfamily\footnotesize,breaklines=true,columns=fullflexible,keepspaces=true,frame=single,rulecolor=\color{light},backgroundcolor=\color{light},xleftmargin=4pt,xrightmargin=4pt}
\newcolumntype{Y}{>{\raggedright\arraybackslash}X}
\begin{document}
\thispagestyle{empty}
{\color{teal}\large ENGINEERING PROJECT REPORT}\par\vspace{8mm}
{\color{navy}\Huge\bfseries DNA Data Storage\\[4pt]Encoder}\par\vspace{5mm}
{\Large Implementation, evaluation,\\testing and web demo}\par\vspace{5mm}
{\large 9 October 2026}\par
\vspace{6mm}
\begin{tabularx}{\linewidth}{YYY}
\toprule
\textcolor{teal}{\Large\bfseries 480 trials} & \textcolor{teal}{\Large\bfseries @@N_TESTS@@ tests} & \textcolor{teal}{\Large\bfseries @@SITE_TRIALS@@ trials}\\
Matched-budget codec validation & Automated checks, @@COVERAGE@@\% line coverage & Error-tolerant web codec benchmark\\
\bottomrule
\end{tabularx}\par
\vspace{6mm}
\textbf{Executive summary.} This project converts opaque file bytes into indexed DNA
oligonucleotides and reconstructs the file from clean or simulated sequencing reads.
It combines several payload mappings with a common container, index, integrity checks,
consensus decoder and outer error correction. A matched-budget experiment runner now
compares codecs at exactly equal total synthesized nucleotide counts.

Five fresh examples were run before preparing this report. Text and binary files were
recovered exactly; a binary file was also recovered after two oligos were removed and
under a 3\% insertion/deletion/substitution (IDS) model at depth 12. A 9\% IDS example
at depth 2 failed with an incomplete header and returned no recovered data.

The earlier 480-trial validation and completed 100 KB clustering experiment provide
broader implementation evidence. They do not establish a general advantage for any
codec, and the full publication study remains unfinished.

An adversarial testing campaign then ran @@N_TESTS@@ Python tests and @@NODE_TOTAL@@ browser-codec
tests. Across corrupted reads, mixed pools, truncated archives and misuse, no decoder ever
reported success with wrong bytes, but @@N_BUGS@@ reproducible defects were found and are
documented as failing tests. A public web demo
(\url{https://dadabsdk.github.io/dna-encoder/}) adds an error-tolerant browser format. In
@@SITE_TRIALS@@ seeded trials it recovered @@REAL_EXACT@@ of @@REAL_TRIALS@@ files exactly at 1\% per-base
error with three reads per strand and 25\% parity, and returned partial files with flagged
uncertain bytes when damage exceeded the parity budget.

\note{\textbf{Evidence boundary.} Fresh examples are single realizations, not estimates of
recovery probability. Historical experiments are identified separately. All successful
fresh examples require both the decoder's integrity checks and exact input/output byte
and SHA-256 equality. No physical DNA synthesis experiment was performed here.}
\vfill
{\small\textbf{Contents:} architecture; physical format and metrics; worked examples;
matched-budget validation; clustering at scale; earlier research phases; testing campaign;
web demo; limitations and next steps; reproduction and evidence.}

\pageheading{Project objective and architecture}{Question: do sequence constraints justify the nucleotides they consume?}
The central design question is whether spending nucleotides on sequence constraints
improves file recovery more than spending the same budget on additional error correction.
The implementation is an integration and trade-off study. It does not claim to invent
rotating ternary mapping, steering redundancy, Reed--Solomon coding or LT fountain coding.

\subsection*{End-to-end workflow}
\begin{center}
\begin{tikzpicture}[node distance=7mm,box/.style={draw=teal,fill=light,rounded corners=2pt,text width=13.3cm,align=center,minimum height=9mm,font=\small},arr/.style={-{Latex[length=2mm]},thick,navy}]
\node[box] (a) {File bytes $\longrightarrow$ type detection and conditional Zstandard compression};
\node[box,below=of a] (b) {Outer RS or LT redundancy $\longrightarrow$ index-bound CRC-32 $\longrightarrow$ payload mapping};
\node[box,below=of b] (c) {Primers + shared index/seed + encoded payload $\longrightarrow$ fixed-length oligo pool};
\node[box,below=of c] (d) {Read channel $\longrightarrow$ orientation and primer trimming $\longrightarrow$ clustering/consensus};
\node[box,below=of d] (e) {Payload decoding and optional repair $\longrightarrow$ outer decoding $\longrightarrow$ file integrity};
\draw[arr] (a)--(b);\draw[arr] (b)--(c);\draw[arr] (c)--(d);\draw[arr] (d)--(e);
\end{tikzpicture}
\end{center}

\subsection*{Payload mappings}
{\small\begin{tabularx}{\linewidth}{p{31mm}Y}
\toprule Mapping & Role and implementation\\\midrule
Naive 2-bit & Quaternary mapping; optional whitening creates the unconstrained whiten+RS comparison.\\
Goldman & Rotating ternary mapping. Uses this project's common ECC, not the original historical overlap scheme.\\
Steering A & Quaternary information symbols with inserted steering positions and seed retries.\\
Steering B & Rotating ternary information symbols with steering; lower raw information rate.\\
Steering C & Run-length-limited quaternary information mapping with steering and a capacity margin.\\
Fountain & Screened LT droplets replace outer RS. CRC and file integrity checks remain in place.\\
\bottomrule\end{tabularx}}

\subsection*{Code organization}
\texttt{encoder.py} and \texttt{decoder.py} orchestrate the pipeline. The
\texttt{codecs/} package contains mappings; \texttt{ecc.py} implements CRC/RS;
\texttt{container.py} handles compression and metadata. \texttt{channel.py} provides
controlled IDS simulations; \texttt{nanopore.py} integrates external channels.
\texttt{benchmark.py} performs budget matching and repeatable evaluation.

\pageheading{Physical format, protection and metrics}{All comparisons must account for the complete physical pool.}
\subsection*{Default oligo layout}
\begin{center}
\begin{tabular}{ccccc}
\toprule Forward primer & Index & Seed & Payload region & Reverse-primer tail\\
20 nt & 15 nt & 3 nt & 142 nt & 20 nt\\\bottomrule
\end{tabular}
\end{center}
The total is 200 nt. The physical tail is the reverse complement of the reverse primer.
The codec-independent index and seed use rotating ternary encoding. Header oligos use
reserved indices below 64 and a fixed bootstrap format, allowing recovery of codec
parameters before data-payload decoding. The default uses five distinct copies of each
header fragment; the number of fragments depends on codec metadata.

\subsection*{Integrity and recovery}
Each data frame has a 32-bit CRC bound to its index. In the default geometry, whitened
2-bit mapping carries 252 data bits plus a 32-bit CRC in 142 nt; Goldman carries 193
data bits plus the same CRC. The CRC's nucleotide cost depends on the mapping.

For RS arms, bytes form rows of width $\lceil d/8\rceil$, with $d$ payload bits per
oligo. RS columns are coded over GF(256), with at most 255 symbols per block. Rows are
bit-packed into oligos, so a missing physical oligo can cross row boundaries. In the
fresh examples, the parity rule is $p=\lceil0.15k\rceil$; parity rows and physical parity
oligos are therefore not interchangeable counts.

\begin{description}
\item[Erasure strength.] CRC-failing payloads are discarded; the outer RS path uses
  erasures and its consistency checks.
\item[Repair strength.] Adds bounded, CRC-guided homopolymer run-length edits and
  outer RS errors-and-erasures decoding. It is the CLI default. Both strengths are
  evaluated on identical read sets in this report.
\end{description}
Decompression is followed by original-length and file CRC-32 verification. For the
worked examples, the caller additionally checks byte equality and SHA-256. These checks
are software integrity measures, not cryptographic authentication of a DNA pool.

\subsection*{Density definitions}
For $N$ fixed-length oligos of length $L$, the charged synthesis budget is
\[
B_{\mathrm{nt}}=\sum_{i=1}^{N} L_i=NL.
\]
The whole-pool densities reported by this implementation are
\[
D_{\mathrm{file}}=\frac{8\,n_{\mathrm{original}}}{B_{\mathrm{nt}}},\qquad
D_{\mathrm{stored}}=\frac{8\,n_{\mathrm{stored}}}{B_{\mathrm{nt}}}.
\]
The separate code-density metric divides stored bits by the data rows' share of
payload-region nucleotides; that denominator includes frame CRC, row padding and steering.
Its accounting is defined in \texttt{dnastore/metrics.py}. Tiny files are dominated by
header overhead. Compression can improve file density without improving the DNA mapping.

\pageheading{Fresh worked examples: recovery outcomes}{Executed for this report using scripts/report\_examples.py.}
\subsection*{Inputs and encoder settings}
The text example contains 205 UTF-8 bytes, including accented Latin and Devanagari text.
Goldman encoding compresses it to 180 stored bytes. The binary input is 1,024 deterministic
pseudorandom bytes generated with NumPy seed 20261008 and uses whitened 2-bit mapping.
Compression is not selected for that input. Both pools use global encoder seed 42,
200-nt oligos, default primers, parity ratio 150 per thousand and one worker.

{\small\begin{tabularx}{\linewidth}{Yrr}
\toprule Quantity & Text pool & Binary pool\\\midrule
Original / stored bytes & 205 / 180 & 1,024 / 1,024\\
Header / data-plus-parity oligos & 10 / 11 & 15 / 38\\
Total oligos / nucleotides & 21 / 4,200 & 53 / 10,600\\
RS data rows / parity rows & 8 / 2 & 32 / 5\\
Whole-pool file density (bits/nt) & 0.3905 & 0.7728\\
Oligos with hard-constraint violations & 0 & 35\\\bottomrule
\end{tabularx}}

The binary pool is deliberately an unconstrained baseline. Its 35 violating oligos
are recorded, not hidden. The examples use the research API; the CLI would require
\texttt{--allow-violations} to write that pool. This does not weaken CRC or file checks.

\subsection*{Observed outcomes}
{\small\begin{tabularx}{\linewidth}{lYrrrr}
\toprule ID & Read condition & Bytes & Reads & Erasure & Repair\\\midrule
@@EXAMPLE_ROWS@@
\bottomrule\end{tabularx}}\par
\smallskip
\textbf{Pass} means exact reconstruction of the original bytes. Both strengths use the
same reads in each row. In E3, indices 64 and 65 were removed; all header oligos were
retained. Thus E3 demonstrates recovery of two data-oligo losses, not arbitrary loss of
two oligos anywhere in the pool.

For E4 and E5, coverage is lognormal with $\sigma=0.27$, orientation is reversed with
probability 0.5, and there is no explicit dropout term. Read allocation is multinomial
at a fixed total count. E4 has $p_s=p_i=p_d=0.01$, depth 12, channel seed 101. E5 has
$p_s=p_i=p_d=0.03$, depth 2, channel seed 202.

\note{\textbf{Failure is part of the demonstration.} E5 reports \texttt{header incomplete}.
The decoder returns no reconstructed file, so data decoding is not established for that
case. Changing both depth and error rate does not identify a recovery threshold or
isolate either factor's causal effect.}

\pageheading{Fresh examples: error burden and byte-level evidence}{Single realizations; descriptive measurements without recovery-probability claims.}
\begin{center}
\begin{tikzpicture}
\begin{axis}[width=.92\linewidth,height=58mm,ybar stacked,bar width=22pt,
 symbolic x coords={E4,E5},xtick=data,ymin=0,ymax=10,
 ylabel={Simulated events per template nt (\%)},
 xticklabels={E4: 3\% IDS / depth 12,E5: 9\% IDS / depth 2},
 tick label style={font=\small},label style={font=\small},
 legend style={at={(.5,1.03)},anchor=south,legend columns=3,font=\small,draw=none},
 ymajorgrids=true,grid style={gray!20}]
\addplot[fill=teal,draw=teal] coordinates {@@SUB_COORDS@@};
\addplot[fill=navy!75,draw=navy!75] coordinates {@@INS_COORDS@@};
\addplot[fill=orange!70,draw=orange!70] coordinates {@@DEL_COORDS@@};
\legend{Substitutions,Insertions,Deletions}
\end{axis}\end{tikzpicture}
\end{center}
{\small\begin{tabularx}{\linewidth}{Yrr}
\toprule Measured quantity & E4 & E5\\\midrule
Substitution events & 1,240 & 620\\
Insertion events & 1,226 & 623\\
Deletion events & 1,257 & 647\\
Total events / original template bases & 3,723 / 127,200 & 1,890 / 21,200\\
Empirical event fraction & 2.927\% & 8.915\%\\\bottomrule
\end{tabularx}}
The denominator is read count times 200 original template bases, including primers.
These simulator event counts are not alignment-derived edit distances or consensus
error rates. E4 has more absolute errors because it has six times as many reads. Both
strengths recover all 38 data-plus-parity oligos in E4; neither requires a successful
CRC-guided payload repair for this particular realization.

\subsection*{One actual encoded oligo}
The first data oligo of the text pool (index 64), wrapped at 50 bases:
\begin{lstlisting}
@@DNA@@
\end{lstlisting}
\subsection*{Input/output SHA-256 checks}
For successful cases, the respective input and output digests are identical. Line
breaks below are for display only; each digest contains 64 hexadecimal characters.
\begin{lstlisting}
@@HASHES@@
\end{lstlisting}

\pageheading{Recorded matched-budget validation}{Historical run: validation\_optimized\_20261008; not rerun for this report.}
This run uses two independent synthetic payload seeds, each 800 bytes, and six arms.
All arms have exactly 62 oligos of 200 nt: \textbf{12,400 total nt per arm}.
Each payload/arm receives 20 independent IDS channel realizations at mean depth 8,
with $p_s=p_i=p_d=0.01$ and lognormal coverage $\sigma=0.27$. Each read set is decoded at
both strengths: $2\times6\times20\times2=480$ decoder trials, corresponding to 240 read sets.
All 12 encoded pools pass noiseless reconstruction. [P3]

{\small\begin{tabularx}{\linewidth}{Yrrrr}
\toprule & \multicolumn{2}{c}{Payload seed 0} & \multicolumn{2}{c}{Payload seed 1}\\
\cmidrule(lr){2-3}\cmidrule(lr){4-5}
Arm & Erasure & Repair & Erasure & Repair\\\midrule
@@VALIDATION_ROWS@@
\bottomrule\end{tabularx}}

\subsection*{Uncertainty and interpretation}
The 95\% Wilson intervals for the counts appearing above are:
\begin{center}\begin{tabular}{lrr}
\toprule Success count & Estimated recovery & Wilson 95\% interval\\\midrule
20/20 & 100\% & 83.89--100.00\%\\
19/20 & 95\% & 76.39--99.11\%\\
17/20 & 85\% & 63.96--94.76\%\\\bottomrule
\end{tabular}\end{center}
For $k$ successes out of $n$ trials and $z=1.96$, the interval is
\[
\frac{\hat p+z^2/(2n)\;\pm\;z\sqrt{\hat p(1-\hat p)/n+z^2/(4n^2)}}{1+z^2/n},
\qquad\hat p=k/n.
\]
Twenty successes do not prove a 100\% population recovery probability. Intervals are
conditional on each fixed payload and configuration; shared-payload trials are not
pooled into a single confidence interval across payloads.

Repair increases steering-B recovery from 17/20 to 20/20 on payload 0, but does not
improve its 19/20 outcome on payload 1. The other five arms succeed in every tested
condition here. This is evidence of correct operation and a finite-sample difference
at one depth, not a general ranking of coding strategies.

\note{\textbf{Budget matching matters.} Extra parity is assigned to denser mappings until
all physical pool sizes match. Header differences, padding and primers are included;
no unused filler oligos are added. Global matching lets the lowest-rate arm influence
the common budget. Pairwise steering-versus-whiten+RS experiments are still needed
across periods, weights, payload types and channel settings.}

\pageheading{Recorded 100 KB clustering-scale experiment}{Historical Phase 3 run; single payload and unmatched synthesis budgets.}
The completed run encodes 100,000 bytes separately with six arms, producing thousands
of oligos per pool. It measures all distinct-body pairs, clusters simulated reads and
attempts complete file recovery at requested mean depth 10. Both calibrated Badread
and a controlled IDS9 channel are evaluated. [P4]

{\small\begin{tabularx}{\linewidth}{Yrrr}
\toprule Arm & Oligos & Min. body ED/$L$ & IDS9 split (\%)\\\midrule
@@SCALE_ROWS@@
\bottomrule\end{tabularx}}
$L=160$ nt is the body length after removing primers. The minimum distances are
measured over all pairs of distinct oligo bodies. No pair is within the $0.30L$
read-clustering radius. This finite-pool observation is not a proof about arbitrary
future pools or error realizations.

\subsection*{Read grouping and file recovery}
\begin{itemize}
\item All 12 arm/channel rows report zero observed cluster merges and zero
  \texttt{groups\_multi\_index} after decoding.
\item Badread has zero measured split rate; all six arms recover at both decoder strengths.
\item IDS9 splits between 2.04\% and 6.41\% of observed oligos across multiple
  multi-read groups. No arm recovers the entire file at either strength.
\end{itemize}
Here IDS9 uses \textbf{4\% substitutions, 2.5\% insertions and 2.5\% deletions}.
It differs from the equal-thirds 9\% IDS model used in fresh example E5. The scale run
uses lognormal coverage spread 0.27 and requested mean depth 10. Pool sizes and synthesis
budgets are not equal across arms; it is a clustering check, not a fair density-versus-
recovery comparison.

\subsection*{Channel scope}
The current project design uses identity-matched Badread as the main read-level
channel and uniform IDS as a controlled reference. The raw-signal
squigulator--Dorado path is retained as a stress test because the project's calibration
found it miscalibrated against available real R10.4.1 storage data. [P5]

The real calibration source consists of plasmid fragments rather than a synthesized
200-nt oligo pool. Badread does not reproduce the small systematic inter-read component
seen in those data. Consequently, simulated recovery is not a wet-lab performance
guarantee. This report summarizes the recorded project calibration; it does not
independently revalidate the external dataset or literature.

\pageheading{Earlier research phases: constraints, structure and channel}{Phase 2--3 measurements. The density sweep predates the format-v3 inner CRC-32 and is not directly comparable with v3 results.}
\subsection*{What sequence constraints cost (Phase 2 sweep, format v2)}
Each arm encoded the same input into 200-nt oligos. Code density counts payload bits per
payload nucleotide; effective density counts file bits per synthesized nucleotide,
including primers, index, header oligos and parity. A violating oligo breaks at least one
hard synthesis constraint (homopolymer, GC content or primer match).
\begin{center}\small
\begin{tabular}{lrrrr}
\toprule
Arm & Code density (b/nt) & Effective density (b/nt) & Max run & Violating oligos\\
\midrule
@@P2_ROWS@@
\bottomrule
\end{tabular}
\end{center}
Unconstrained quaternary mappings are densest but most of their oligos break hard
constraints; whitening alone does not fix this. Every constrained arm met the hard
constraints at a measurable density cost. Whether that cost buys more recovery than the
same nucleotides spent on parity is the open question of the study.

\subsection*{A hidden structure cost of homopolymer limits}
Median minimum free energy (kcal/mol, 95\% CI) of 2{,}000 random 160-nt bodies; more
negative means more stable self-folding.
\begin{center}\small
\begin{tabular}{lrrr}
\toprule
Sequence model & MFE, 37\,\textdegree C & MFE, 60\,\textdegree C & Paired fraction\\
\midrule
@@MFE_ROWS@@
\bottomrule
\end{tabular}
\end{center}
Run-length-1 coding folds markedly more stably. Banning one successor raises the
probability of reverse-complementary k-mer pairs, and a control that bans a different
successor reproduces most of the gap. Constrained alphabets therefore start from a
structural deficit that their steering must repay.

\subsection*{Choosing a structure metric}
Spearman correlation with equilibrium primer-site accessibility at 60\,\textdegree C on a
held-out half (n = 505):
\begin{center}\small
\begin{tabular}{lr}
\toprule
Proxy & $\rho$\\
\midrule
@@PROXY_ROWS@@
\bottomrule
\end{tabular}
\end{center}
Global MFE is a poor surrogate for primer accessibility, so accessibility became the
primary structure metric; the cheap k-mer proxy is only a tie-breaker.

\subsection*{Read-channel calibration against real nanopore data}
Real reads: public R10.4.1 data (Chen et al.\ 2025, Zenodo 10.5281/zenodo.16883332),
Dorado basecalls on identical reference intervals.
\begin{center}\small
\begin{tabular}{lrrrr}
\toprule
Measure & Real hac & Real sup & Squigulator $\to$ Dorado & Badread (matched)\\
\midrule
@@CAL_ROWS@@
\bottomrule
\end{tabular}
\end{center}
Identity-matched Badread reproduces the real error level and is the main read-level
channel; squigulator is retained only as a labelled, miscalibrated stress channel.
Badread omits the systematic read-to-read errors seen in real data.

\pageheading{Testing campaign and known defects}{Full suites rerun for this report; evidence in results/testing\_20261009/.}
\begin{center}\small
\begin{tabular}{lr}
\toprule
Python tests (pytest) & @@N_TESTS@@\\
\quad passed / skipped / failed & @@N_PASS@@ / @@N_SKIP@@ / @@N_FAIL@@\\
\quad known defects (strict expected failures) & @@N_XFAIL@@\\
Browser-codec tests (Node.js) passed / failed & @@NODE_PASS@@ / @@NODE_FAIL@@\\
Line coverage of \texttt{dnastore/} & @@COVERAGE@@\%\\
\bottomrule
\end{tabular}
\end{center}
\textbf{What was attacked.} Heavy insertion/deletion/substitution noise, truncated and
chimeric reads, junk-only input, pools mixing two files with the same primers, outer
Reed--Solomon erasures at and beyond capacity (checked against the \texttt{reedsolo}
reference), undetected corrupted oligos, 556 corrupted or truncated stream archives,
33 malformed FASTA/FASTQ cases, determinism across hash seeds and worker counts, and the
local web server's host, token and path-traversal checks. @@N_GUARDS@@ guard tests now
pin these properties. \textbf{No test produced a successful decode with wrong bytes.}

\textbf{Defects found.} Each is a strict expected-failure test in
\texttt{tests/test\_adversarial.py}: the suite stays green, and a fix makes the test pass,
which then fails the run until the marker is removed. None is fixed yet.
\begin{center}\small
\begin{tabularx}{\linewidth}{lY}
\toprule
ID & Defect\\
\midrule
@@BUG_ROWS@@
\bottomrule
\end{tabularx}
\end{center}
\note{\textbf{Impact.} The fountain-codec defect affects benchmark interpretation for small
files: a pool that cannot be decoded from perfect reads makes later channel trials
meaningless. The benchmark records \texttt{noiseless\_ok} per pool, which must be checked.
The \texttt{reedsolo} defect breaks the default CLI decoder on a plain install.}

\pageheading{Web demo and the error-tolerant browser codec}{\url{https://dadabsdk.github.io/dna-encoder/}; benchmark: scripts/site\_resilience.mjs.}
The static site runs entirely in the browser and is deployed by GitHub Actions with
commit-versioned asset URLs. Its default format (\texttt{DNASTORE-RS-1}) is separate from
the research codec: 32-byte payloads become 160-nt strands
[index $\|$ payload $\|$ CRC-32] with keystream whitening, plus Reed--Solomon parity
strands across each block of up to 255 strands and eight copies of the metadata. The
decoder (i) accepts reads whose CRC verifies, (ii) repairs one substitution, insertion or
deletion per read by CRC-guided search, (iii) merges several reads of one strand by
banded-alignment consensus, (iv) rebuilds missing strands from parity, and (v) beyond the
parity budget returns the file with unrecoverable bytes reported as uncertain. Bytes not
reported uncertain are exact; the benchmark asserts this in every trial.

\textbf{Benchmark.} @@SITE_TRIALS@@ trials: random 20{,}000-byte payloads, @@SITE_SEEDS@@ seeds per
cell; 5\% strand dropout; per-base error split 50\% substitution, 25\% insertion, 25\%
deletion, independent per read. Synthesis cost: @@NTB_100@@ / @@NTB_250@@ / @@NTB_500@@ nt per byte at
10 / 25 / 50\% parity (raw 2-bit mapping: 4 nt per byte).
\begin{center}
\begin{tikzpicture}
\begin{axis}[width=0.62\linewidth,height=5.2cm,xlabel={Reads per strand},ylabel={Exact recoveries (of @@SITE_SEEDS@@)},
  xtick={1,2,3},ymin=-0.5,ymax=10.5,legend pos=outer north east,legend style={font=\footnotesize},title={25\% parity},
  title style={font=\small},grid=major,grid style={gray!25}]
@@SITE_PLOT@@
\end{axis}
\end{tikzpicture}
\end{center}
\par\noindent\begin{minipage}{\linewidth}
\textbf{Exact recoveries by parity level.}
\begin{center}\small
\begin{tabular}{rrrrr}
\toprule
Error/base & Reads & 10\% parity & 25\% parity & 50\% parity\\
\midrule
@@SITE_GRID@@
\bottomrule
\end{tabular}
\end{center}
\end{minipage}
\par\noindent\begin{minipage}{\linewidth}
\textbf{Partial recovery at 25\% parity.}
\begin{center}\small
\begin{tabular}{rrrrrr}
\toprule
Err. & Reads & Exact & Guaranteed & Correct & ms\\
\midrule
@@SITE_PARTIAL@@
\bottomrule
\end{tabular}
\end{center}
\end{minipage}
{\small \emph{Guaranteed} is the mean share of bytes not flagged
uncertain; \emph{Correct} also counts flagged bytes whose best guess was right; ms is the
median decode time (Node.js @@NODE_VERSION@@).}

Multiple reads per strand matter more than parity: a single read at 1--2\% error per base
usually carries two or more errors, which single-edit repair cannot fix. The channel is
synthetic with independent errors; real systematic errors and PCR bias would lower these
numbers. This browser format is a teaching and robustness demonstration, not a replacement
for the constrained research codec, and is not yet readable by the CLI.

\pageheading{Engineering quality and remaining limitations}{The software workflow is operational; the broader research study is not complete.}
\subsection*{Completed reliability work}
\begin{itemize}
\item Added the missing matched-budget benchmark command, independent channel seeds,
  both decoder strengths, Wilson intervals, empirical depth search, output tables,
  plots, source snapshots and provenance records.
\item Added multiline FASTA, wrapped FASTQ and gzip support with malformed-record checks.
  Ambiguous reads containing \texttt{N} are discarded and counted; lowercase ACGT is normalized.
\item Validated key RS, header, geometry, probability and codec parameters. Made the
  strict CLI constraint policy apply to baseline codecs as well.
\item Bundled default primers so installed-wheel roundtrips work outside the checkout.
\item Excluded groups composed entirely of verified header candidates from futile
  data-repair attempts, while preserving uncertain and mixed header/data clusters.
\end{itemize}
The prior expanded full suite passed \textbf{240 tests with one skip}. After the header
optimization, the affected suites passed 44 tests; header-only and mixed-cluster
regressions also passed separately. These are overlapping verification runs, not
additive test counts. Dependency consistency, wheel construction and an installed-wheel
roundtrip were also checked. [P1]

A saved single-case profile showed repair attempts falling from 17 to 5 and profiled
time from 9.32 to 3.53 seconds. Across 120 shared before/after trials, recovery and
continuous error metrics were identical. Timing was measured with concurrent background
load; it is not a general throughput claim.

\subsection*{Limitations}
\begin{enumerate}
\item \textbf{Evidence size.} Fresh examples use one seed per noisy setting; the
  matched validation uses only two small synthetic payloads at one depth.
\item \textbf{Physical fidelity.} IDS is a controlled model. Badread lacks systematic
  inter-read errors. Neither proves synthesis, PCR or long-term storage reliability.
\item \textbf{Constraint costs.} Baselines may emit violating sequences in research
  mode. The synthesis-side penalty of these violations is not established by this report.
\item \textbf{Inference.} The bisection depth result is an empirical boundary. Independent
  stochastic outcomes need not be monotone; right censoring and observed reversals must be reported.
\item \textbf{Operations.} Output directories are protected from overwrite. Complete
  depth probes are saved incrementally; automatic interrupted-run resume is not implemented.
\end{enumerate}
\textbf{Next steps.} Run the larger \texttt{configs/paper.yaml} study, add a licensed
real-file corpus, sweep steering periods and normalized weights, and perform the specified
coverage/accessibility sensitivities with explicit assumption labels. Complete literature
verification before making novelty or publication-level claims.

\pageheading{Reproduction and evidence manifest}{All paths below are relative to the project root unless stated otherwise.}
\subsection*{Reproduce the fresh demonstrations}
\begin{lstlisting}
.venv/bin/python scripts/report_examples.py \
  --out results/report_examples_rerun
\end{lstlisting}
The output directory must be new or empty. The script saves original files, physical
pools, exact read sets, successful reconstructions, simulator counts, decoder reports,
configuration, SHA-256 digests and source checksums. E5 intentionally has no recovered file.

\subsection*{Run the benchmark and compile this document}
\begin{lstlisting}
.venv/bin/dnastore benchmark configs/benchmark_validation.yaml \
  --out results/phase4/validation_rerun
.venv/bin/dnastore benchmark configs/paper.yaml \
  --out results/phase4/full_study
.venv/bin/python scripts/build_project_report.py
pdflatex -interaction=nonstopmode -halt-on-error \
  -output-directory=output/pdf output/report/dna_storage_project_report.tex
pdflatex -interaction=nonstopmode -halt-on-error \
  -output-directory=output/pdf output/report/dna_storage_project_report.tex
cp output/pdf/dna_storage_project_report.pdf reports/
\end{lstlisting}
This source is standalone: diagrams and chart coordinates are embedded; no external
images or bibliography files are required. A standard LaTeX installation with TikZ,
PGFPlots, Latin Modern, microtype and the listed table/layout packages is sufficient.
The full study can be substantially more expensive than the demonstrations.

\subsection*{Rerun the tests and the web-codec benchmark}
\begin{lstlisting}
.venv/bin/pytest -q --cov=dnastore --cov-report=json:OUT/coverage.json \
  --junitxml=OUT/pytest_junit.xml
node --test --test-reporter=tap site/tests/resilient.test.mjs
node scripts/site_resilience.mjs results/site_resilience_rerun
\end{lstlisting}

\subsection*{Primary project evidence}
{\small\begin{description}
\item[P1] \path{docs/PROJECT_STATUS.md} and \path{docs/DESIGN.md}: design decisions,
  completed engineering work and recorded verification history.
\item[P2] \path{results/report_examples_20261008/examples.json} and \path{results.csv}
  in the same directory: fresh examples, complete digests and decoder reports.
\item[P3] \path{results/phase4/validation_optimized_20261008/}: \path{summary.csv},
  \path{trials.csv}, pool manifests and \path{provenance.json}; historical 480-trial validation.
\item[P4] \path{results/phase3/cluster_scale_100KB/cluster_scale.csv}: all 12
  historical scale-run arm/channel rows.
\item[P5] \path{docs/notes/channel_calibration.md}: calibration methods, recorded
  channel measurements and transfer limitations.
\item[P6] \path{docs/BENCHMARK.md}, \path{dnastore/metrics.py} and
  \path{dnastore/benchmark.py}: budget, metric and statistical definitions.
\item[P7] \path{results/phase2/sweep_v2/summary.csv}, \path{results/phase2/composition/mfe_summary.csv}
  and \path{docs/notes/structure.md}: historical constraint, structure and proxy results.
\item[P8] \path{results/testing_20261009/}: JUnit XML, coverage JSON and Node TAP output of
  the test rerun; defects parsed from \path{tests/test_adversarial.py}.
\item[P9] \path{results/site_resilience_20261009/}: per-trial \path{results.csv},
  \path{summary.csv} and \path{provenance.json} (codec SHA-256 recorded).
\end{description}}
\note{\textbf{Traceability.} The report builder reads the saved JSON/CSV evidence rather
than rerunning or fabricating measurements. \texttt{output/report/report\_evidence.json}
records SHA-256 hashes of the evidence files and generated LaTeX source. Rebuilding with
changed evidence requires reviewing the interpretation as well as the numbers.}
\end{document}
'''
for key,value in {'EXAMPLE_ROWS':'\n'.join(rows),'VALIDATION_ROWS':'\n'.join(vt),'SCALE_ROWS':'\n'.join(st),
                  'DNA':'\n'.join(dna[i:i+50] for i in range(0,len(dna),50)),
                  'HASHES':hashes,'SUB_COORDS':errors[0],'INS_COORDS':errors[1],'DEL_COORDS':errors[2],
                  'P2_ROWS':'\n'.join(p2),'MFE_ROWS':'\n'.join(mfe_rows),'PROXY_ROWS':'\n'.join(proxy_rows),'CAL_ROWS':'\n'.join(cal_rows),
                  'N_TESTS':str(len(cases)),'N_PASS':str(n_pass),'N_SKIP':str(n_skip),'N_FAIL':str(n_fail),'N_XFAIL':str(n_xfail),
                  'NODE_PASS':str(node_pass),'NODE_FAIL':str(node_fail),'NODE_TOTAL':str(node_pass+node_fail),'COVERAGE':f'{coverage:.0f}',
                  'N_BUGS':str(len(bugs)),'N_GUARDS':str(guards),'BUG_ROWS':'\n'.join(bug_rows),
                  'SITE_TRIALS':f'{site_trials:,}','SITE_SEEDS':site[0]['trials'],'SITE_PLOT':site_plot,'SITE_GRID':'\n'.join(site_grid),
                  'SITE_PARTIAL':'\n'.join(site_partial),'NODE_VERSION':esc(site_prov['node']),
                  'NTB_100':f"{nt_per_byte[100]:.2f}",'NTB_250':f"{nt_per_byte[250]:.2f}",'NTB_500':f"{nt_per_byte[500]:.2f}",
                  'REAL_EXACT':real['exact'],'REAL_TRIALS':real['trials']}.items():
    tex=tex.replace('@@'+key+'@@',value)
assert '@@' not in tex
out=ROOT/'output/report/dna_storage_project_report.tex';out.parent.mkdir(parents=True,exist_ok=True);out.write_text(tex)
files=[EX/'examples.json',EX/'results.csv',VALID/'summary.csv',VALID/'provenance.json',ROOT/'results/phase3/cluster_scale_100KB/cluster_scale.csv',
       SWEEP,MFE,STRUCT,CAL,TEST/'pytest_junit.xml',TEST/'coverage.json',TEST/'node_tests.tap',ADV,
       SITE/'results.csv',SITE/'summary.csv',SITE/'provenance.json',ROOT/'site/resilient.js',
       ROOT/'docs/PROJECT_STATUS.md',ROOT/'docs/DESIGN.md',out]
manifest={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
(out.parent/'report_evidence.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(out)
