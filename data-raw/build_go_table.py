import gzip, sys, os, collections, re
# Usage: python3 data-raw/build_go_table.py go-basic.obo HUMAN-uniprot.gaf.gz OUTDIR [--app-dir data/go_release]
#   OUTDIR receives the full tables (all_terms, term2gene_long, BP/MF/CC .gmt).
#   --app-dir also writes the three files the app's "GO (official release)" databases read.
obo, gaf, outdir = sys.argv[1:4]
app_dir = sys.argv[sys.argv.index("--app-dir") + 1] if "--app-dir" in sys.argv else None
EXP = {"EXP","IDA","IPI","IMP","IGI","IEP","HTP","HDA","HMP","HGI","HEP"}
# ---- parse ontology
terms = {}; alt = {}; parents = collections.defaultdict(set); cur = None; dv = None
with open(obo) as f:
    for line in f:
        line = line.rstrip("\n")
        if dv is None and line.startswith("data-version:"): dv = line.split(":",1)[1].strip()
        if line.startswith("["):
            cur = {} if line == "[Term]" else None
            continue
        if cur is None or not line: continue
        k, _, v = line.partition(": ")
        if k == "id": cur["id"] = v; terms[v] = cur; cur.update(obsolete=False)
        elif k == "name": cur["name"] = v
        elif k == "namespace": cur["ns"] = v
        elif k == "is_obsolete" and v == "true": cur["obsolete"] = True
        elif k == "alt_id": alt[v] = cur["id"]
        elif k == "is_a": parents[cur["id"]].add(v.split()[0])
        elif k == "relationship":
            rel, tgt = v.split()[:2]
            if rel == "part_of": parents[cur["id"]].add(tgt)
live = {t for t, d in terms.items() if t.startswith("GO:") and not d["obsolete"]}
anc_cache = {}
def ancestors(t):
    if t in anc_cache: return anc_cache[t]
    out = {t}
    for p in parents.get(t, ()):
        if p in live: out |= ancestors(p)
    anc_cache[t] = out; return out
sys.setrecursionlimit(10000)
# ---- parse annotations
direct = collections.defaultdict(set); direct_exp = collections.defaultdict(set)
hdr = {}; n_lines = n_not = n_obs = 0
with gzip.open(gaf, "rt") as f:
    for line in f:
        if line.startswith("!"):
            m = re.match(r"!\s*([\w-]+):\s*(.*)", line.strip())
            if m: hdr.setdefault(m.group(1), m.group(2))
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 15 or c[0] != "UniProtKB" or c[12] != "taxon:9606": continue
        n_lines += 1
        if "NOT" in c[3].split("|"): n_not += 1; continue
        go = alt.get(c[4], c[4])
        if go not in live: n_obs += 1; continue
        direct[go].add(c[2])
        if c[6] in EXP: direct_exp[go].add(c[2])
allg = collections.defaultdict(set); allexp = collections.defaultdict(set)
for go, gs in direct.items():
    for a in ancestors(go): allg[a] |= gs
for go, gs in direct_exp.items():
    for a in ancestors(go): allexp[a] |= gs
rel = (dv or "").split("/")[-1]
tag = f"GO_human_{rel}"
nsab = {"biological_process":"BP","molecular_function":"MF","cellular_component":"CC"}
out = os.path.join(outdir, f"{tag}_all_terms.tsv.gz")
with gzip.open(out, "wt") as o:
    o.write("go_id\tname\tontology\tn_genes_all\tn_genes_direct\tn_genes_experimental\tgenes_all\tgenes_direct\tgenes_experimental\n")
    for t in sorted(live):
        d = terms[t]
        ga, gd, ge = sorted(allg.get(t,())), sorted(direct.get(t,())), sorted(allexp.get(t,()))
        o.write("\t".join([t, d.get("name",""), nsab.get(d.get("ns"),d.get("ns","")), str(len(ga)), str(len(gd)), str(len(ge)),
                           ",".join(ga), ",".join(gd), ",".join(ge)]) + "\n")
# long format
with gzip.open(os.path.join(outdir, f"{tag}_term2gene_long.tsv.gz"), "wt") as o:
    o.write("go_id\tontology\tgene_symbol\tdirect\texperimental\n")
    for t in sorted(allg):
        ns = nsab.get(terms[t].get("ns"))
        for g in sorted(allg[t]):
            o.write(f"{t}\t{ns}\t{g}\t{int(g in direct.get(t,()))}\t{int(g in allexp.get(t,()))}\n")
# GMT per ontology
for ns, ab in nsab.items():
    with open(os.path.join(outdir, f"{tag}_{ab}.gmt"), "w") as o:
        for t in sorted(allg):
            if terms[t].get("ns") == ns:
                o.write("\t".join([t, terms[t].get("name","")] + sorted(allg[t])) + "\n")
cnt = collections.Counter(nsab.get(terms[t].get("ns")) for t in live)
ann = collections.Counter(nsab.get(terms[t].get("ns")) for t in allg)
print("data-version", dv); print("gaf header", {k:hdr[k] for k in hdr if k in ("gaf-version","date-generated","go-version","generated-by")})
print("lines", n_lines, "NOT dropped", n_not, "obsolete/unknown dropped", n_obs)
print("live terms", dict(cnt), "terms with >=1 gene", dict(ann))
print("genes", len(set().union(*direct.values())))
for t in ("GO:0060337","GO:0034340","GO:0051607"):
    gs = allg.get(t,set()); print(t, terms[t]["name"], len(gs), "IFIT1/2/3:", [g for g in ("IFIT1","IFIT2","IFIT3") if g in gs])
print("GO:0060337 genes:", ",".join(sorted(allg["GO:0060337"])))

if app_dir:
    os.makedirs(app_dir, exist_ok=True)
    with gzip.open(os.path.join(app_dir, "term2gene.tsv.gz"), "wt") as o:
        o.write("go_id\tontology\tgene_symbol\n")
        for t in sorted(allg):
            ns = nsab.get(terms[t].get("ns"))
            for g in sorted(allg[t]): o.write(f"{t}\t{ns}\t{g}\n")
    with gzip.open(os.path.join(app_dir, "term2name.tsv.gz"), "wt") as o:
        o.write("go_id\tname\tontology\n")
        for t in sorted(allg):
            o.write(f"{t}\t{terms[t].get('name','')}\t{nsab.get(terms[t].get('ns'))}\n")
    with open(os.path.join(app_dir, "VERSION.txt"), "w") as o:
        o.write(f"release: {rel}\n"
                f"ontology: go-basic.obo, data-version {dv}\n"
                f"annotations: {os.path.basename(gaf)}, {hdr.get('generated-by','')}, date-generated {hdr.get('date-generated','')}\n"
                "evidence codes: all (IEA and IBA included); NOT-qualified annotations dropped\n"
                "propagated over: is_a, part_of\n"
                "genes: identified by gene symbol; matched to the loaded DE table by symbol\n"
                "built with: data-raw/build_go_table.py --app-dir data/go_release\n")
    print("app files written to", app_dir)
