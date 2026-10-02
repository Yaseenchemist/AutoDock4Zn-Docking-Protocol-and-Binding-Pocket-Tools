#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
  AUTO-DOCK PIPELINE  (AutoDock Vina)  -- fully automated
--------------------------------------------------------------------------------
  انطي بس الـ PDB ID، والباقي كله اوتوماتيك:
    1) ينزّل الـ PDB من RCSB
    2) يكتشف الكوكرستال ليكاند اوتوماتيك (يرفض الماي/الايونات/الاضافات)
    3) يفصل البروتين وينضّفه + يحضّره PDBQT (prepare_receptor4.py)
    4) يحدد الـ pocket/grid box اوتوماتيك من موقع الكوكرستال
    5) REDOCK للّيكاند الاصلي + يحسب RMSD  (validation للبروتوكول)
    6) يحضّر مركبك من الـ SMILES ويعمل docking
    7) يحلل بطاقات الارتباط (PLIP: H-bonds, hydrophobic, pi-stacking...)
    8) يخزن كلشي بفولدر باسم البروتين:  results/<PDB_ID>/
================================================================================

  الأماكن الوحيدة اللي تعدلها بيها = بلوك CONFIG بس (تحت مباشرة).
  كل شي بعده اوتوماتيك، ما تلمسه.
"""

import os, sys, shutil, subprocess, urllib.request, textwrap
import numpy as np

# ============================================================================ #
# ============================  >>> CONFIG <<<  ============================== #
# ============  (هنا بس تعدل — عادةً تغيّر PDB_ID فقط)  ==================== #
# ============================================================================ #

PDB_ID       = "2WI1"          # <-- (1) بدّل الـ PDB ID مال البروتين كل مرة

QUERY_SMILES = "COC1=C(C(=C2C(=C1)C=CC(=O)O2)OC)O"   # (2) مركبك (ثابت لكل الحالات)

# ---- خيارات اختيارية: خليهن على الديفولت الا اذا احتجت ----
LIGAND_RESNAME = ""       # حدّد اسم الليكاند يدوي (مثل "SO4") لو الأوتوماتيك غلط. فارغ = اوتوماتيك
KEEP_CHAINS    = ""       # خلي chain معين فقط، مثل "A". فارغ = كل السلاسل
BOX_PADDING    = 5.0      # هامش البوكس حول الكوكرستال (Å)
EXHAUSTIVENESS = 32       # دقة البحث (اكبر = ادق واثقل)
NUM_MODES      = 20       # عدد الـ poses المخرجة
SEED           = 42       # للتكرار (reproducibility)
RMSD_CUTOFF    = 2.0      # اذا redock RMSD اصغر منه => البروتوكول موثوق

# ---- مسارات الأدوات: خليهن "" حتى يلكاهن اوتوماتيك من الـ PATH ----
# لو ما لكاهن، احط المسار الكامل يدوي هنا.
VINA_BIN        = ""      # مثال: "/home/USER/SOFTWARE/vina/vina"
MGLTOOLS_PYSH   = "/home/yaseen/SOFTWARE/MGLTools-1.5.6/bin/pythonsh"
PREP_RECEPTOR   = "/home/yaseen/SOFTWARE/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"
PREP_LIGAND     = "/home/yaseen/SOFTWARE/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"

# ============================================================================ #
# ========================  END OF CONFIG — لا تعدل تحت  ===================== #
# ============================================================================ #


# --------------------------- helpers -------------------------------------- #
def log(msg):  print(f"[auto-dock] {msg}", flush=True)
def die(msg):  log("ERROR: " + msg); sys.exit(1)

def run(cmd, cwd=None, quiet=False):
    """Run shell command, raise on failure."""
    if not quiet: log("$ " + (cmd if isinstance(cmd, str) else " ".join(cmd)))
    r = subprocess.run(cmd, cwd=cwd, shell=isinstance(cmd, str),
                       capture_output=True, text=True)
    if r.returncode != 0:
        log("stdout:\n" + r.stdout)
        log("stderr:\n" + r.stderr)
        raise RuntimeError(f"command failed: {cmd}")
    return r.stdout

def find_tool(explicit, names):
    """Return explicit path if set & exists, else search PATH & common dirs."""
    if explicit:
        if os.path.exists(explicit): return explicit
        die(f"ما لكيت الاداة بالمسار: {explicit}")
    for n in names:
        p = shutil.which(n)
        if p: return p
    # common MGLTools install locations
    guesses = []
    home = os.path.expanduser("~")
    for base in (f"{home}/SOFTWARE", f"{home}", "/opt", "/usr/local"):
        for pat in ("mgltools*", "MGLTools*", "ADFRsuite*"):
            import glob
            guesses += glob.glob(f"{base}/{pat}")
    for g in guesses:
        for n in names:
            for sub in ("bin", "MGLToolsPckgs/AutoDockTools/Utilities24"):
                cand = os.path.join(g, sub, n)
                if os.path.exists(cand): return cand
    return None


# --------------------------- tool discovery ------------------------------- #
def discover_tools():
    t = {}
    t["vina"]  = find_tool(VINA_BIN, ["vina", "vina_1.2.5_linux_x86_64", "vina.exe"])
    t["pysh"]  = find_tool(MGLTOOLS_PYSH, ["pythonsh"])
    t["prep_r"]= find_tool(PREP_RECEPTOR, ["prepare_receptor4.py", "prepare_receptor.py"])
    t["prep_l"]= find_tool(PREP_LIGAND,   ["prepare_ligand4.py", "prepare_ligand.py"])
    t["obabel"]= shutil.which("obabel")
    t["obrms"] = shutil.which("obrms")
    t["plip"]  = shutil.which("plip")
    if not t["vina"]:  die("ما لكيت vina — احط مساره بـ VINA_BIN")
    if not (t["prep_r"] and t["pysh"]):
        log("تنبيه: prepare_receptor4.py / pythonsh مو موجودين — راح استخدم obabel fallback")
    if not t["obabel"]: log("تنبيه: obabel مو منصّب — بعض الخطوات راح تفشل. نصّبه: sudo apt install openbabel")
    if not t["obrms"]:  log("تنبيه: obrms مو موجود — راح احسب RMSD بطريقة داخلية")
    if not t["plip"]:   log("تنبيه: plip مو منصّب — تحليل الارتباط راح ينطي. نصّبه: pip install plip")
    log("الأدوات: " + ", ".join(f"{k}={'✓' if v else '✗'}" for k,v in t.items()))
    return t


# --------------------------- PDB download --------------------------------- #
def download_pdb(pdb_id, out_pdb):
    if os.path.exists(out_pdb):
        log(f"موجود مسبقاً: {out_pdb}"); return
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    log(f"تنزيل {url}")
    urllib.request.urlretrieve(url, out_pdb)


# --------------------------- ligand auto-detect --------------------------- #
IGNORE_HET = {
    "HOH","WAT","DOD","H2O",
    "NA","CL","K","MG","CA","ZN","MN","FE","FE2","CU","CO","NI","CD","HG","BR","IOD","Iod","I",
    "SO4","PO4","GOL","EDO","PEG","PG4","PGE","1PE","ACT","DMS","MPD","FMT","CIT","TRS","TAR","TLA",
    "IMD","MES","EPE","BME","NAG","BMA","MAN","FUC","NO3","BCT","SCN","AZI","FLC","OGA","MLI","GAI",
    "PE4","P6G","2PE","ACY","BEN","CAC","DTT","BOG","LDA","MYR","OCT","PLM","STE","UNX","UNK","PEO"
}
COFACTOR_HINT = {"HEM","FAD","FMN","NAD","NAP","ADP","ATP","GTP","GDP","SAM","SAH","PLP","COA","TPP"}

from Bio.PDB import PDBParser, PDBIO, Select

def detect_ligand(pdb_path):
    """Return (chain, resname, resseq) of best co-crystal ligand candidate."""
    s = PDBParser(QUIET=True).get_structure("x", pdb_path)
    cands = {}
    for res in s.get_residues():
        hetflag = res.id[0]
        if hetflag.strip() == "":  continue           # standard residues
        rn = res.resname.strip()
        if rn in IGNORE_HET:       continue
        if KEEP_CHAINS and res.parent.id not in KEEP_CHAINS: continue
        n_heavy = sum(1 for a in res if a.element != "H")
        if n_heavy < 6 and rn not in COFACTOR_HINT:    # too small to be a real ligand
            continue
        cands[(res.parent.id, rn, res.id[1])] = n_heavy
    if LIGAND_RESNAME:
        for k in cands:
            if k[1] == LIGAND_RESNAME: return k, cands
        die(f"ما لكيت الليكاند اللي طلبته: {LIGAND_RESNAME}. المرشحين: {list(cands)}")
    if not cands:
        die("ما گدرت اكتشف كوكرستال ليكاند. حدّد LIGAND_RESNAME يدوي.")
    best = max(cands, key=cands.get)
    return best, cands


# --------------------------- selectors ------------------------------------ #
class ProteinSelect(Select):
    """Keep protein only: drop water/HETATM, strip H, keep first altLoc."""
    def accept_residue(self, res):
        if KEEP_CHAINS and res.parent.id not in KEEP_CHAINS: return 0
        return 1 if res.id[0].strip() == "" else 0
    def accept_atom(self, atom):
        if atom.element == "H": return 0                 # re-added cleanly later
        if atom.get_altloc() not in (" ", "", "A"): return 0
        return 1

class LigandSelect(Select):
    def __init__(self, key): self.key = key
    def accept_residue(self, res):
        return 1 if (res.parent.id, res.resname.strip(), res.id[1]) == self.key else 0

def extract_ligand_pdb(raw_pdb, key, out_pdb):
    """Write the ligand atoms PLUS their CONECT records so bond perception is
    correct (RCSB ships CONECT for ligands; Biopython would drop them, which
    makes obabel fragment big/complex ligands)."""
    chain, resname, resseq = key
    atom_lines, serials = [], set()
    for ln in open(raw_pdb):
        if ln.startswith(("ATOM", "HETATM")):
            if (ln[17:20].strip() == resname and ln[21] == chain
                    and ln[22:26].strip() == str(resseq)):
                serials.add(ln[6:11].strip())
                atom_lines.append(ln)
    conect = []
    for ln in open(raw_pdb):
        if ln.startswith("CONECT"):
            refs = [ln[i:i+5].strip() for i in range(6, len(ln.rstrip("\n")), 5)]
            if refs and refs[0] in serials:
                conect.append(ln)
    with open(out_pdb, "w") as f:
        f.writelines(atom_lines); f.writelines(conect); f.write("END\n")
    return len(atom_lines)

def write_selection(pdb_path, out_path, selector):
    s = PDBParser(QUIET=True).get_structure("x", pdb_path)
    io = PDBIO(); io.set_structure(s); io.save(out_path, selector)


# --------------------------- grid box ------------------------------------- #
def grid_box_from_ligand(ligand_pdb, padding):
    s = PDBParser(QUIET=True).get_structure("l", ligand_pdb)
    coords = np.array([a.coord for a in s.get_atoms() if a.element != "H"])
    center = coords.mean(0)
    size   = (coords.max(0) - coords.min(0)) + 2*padding
    size   = np.maximum(size, 8.0)          # floor so box isn't tiny
    return center, size


# --------------------------- receptor prep -------------------------------- #
def pdbqt_coords_ok(path):
    """True if every ATOM/HETATM line has parseable x/y/z at fixed columns."""
    if not os.path.exists(path): return False
    for ln in open(path):
        if ln.startswith(("ATOM", "HETATM")):
            try:
                float(ln[30:38]); float(ln[38:46]); float(ln[46:54])
            except ValueError:
                return False
    return True

def _rebuild_pdbqt_line(ln):
    """Rebuild a column-broken ATOM/HETATM line. Coords/charge/type are recovered
    from the right (space-separated & reliable); a short safe name is synthesized
    so an over-long atom name can't overflow the columns again."""
    rec = ln[:6]
    body = ln[6:].split()
    try:
        atype  = body[-1]
        charge = float(body[-2]); temp = float(body[-3]); occ = float(body[-4])
        z = float(body[-5]); y = float(body[-6]); x = float(body[-7])
        serial = int(body[0])
    except (ValueError, IndexError):
        return None
    name = atype if len(atype) <= 4 else atype[:4]
    return (f"{rec:<6}{serial:>5} {name:<4} {'UNL':>3} A{1:>4}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{temp:6.2f}    {charge:6.3f} {atype:<2}\n")

def sanitize_pdbqt(path):
    """Reformat any ATOM/HETATM line whose x/y/z columns don't parse (fixes the
    MGLTools long-hydrogen-name overflow). Returns True if anything changed."""
    changed, out = False, []
    for ln in open(path):
        if ln.startswith(("ATOM", "HETATM")):
            try:
                float(ln[30:38]); float(ln[38:46]); float(ln[46:54])
            except ValueError:
                fixed = _rebuild_pdbqt_line(ln)
                changed = True
                if fixed: out.append(fixed)
                continue
        out.append(ln)
    if changed:
        open(path, "w").writelines(out)
    return changed

def pdbqt_ligand_ok(path, min_heavy=None):
    """Ligand PDBQT must parse, have exactly one ROOT, and (if min_heavy given)
    contain enough heavy atoms — rejects partial/fragmented outputs."""
    if not os.path.exists(path): return False
    roots = heavy = 0
    for ln in open(path):
        if ln.startswith("ROOT"): roots += 1
        if ln.startswith(("ATOM", "HETATM")):
            try: float(ln[30:38]); float(ln[38:46]); float(ln[46:54])
            except ValueError: return False
            toks = ln.rstrip().split()
            atype = toks[-1] if toks else ""
            if not atype.upper().startswith("H"): heavy += 1
    if roots != 1: return False
    if min_heavy and heavy < max(3, int(0.7 * min_heavy)): return False
    return True

def prepare_receptor(protein_pdb, out_pdbqt, tools):
    ok = False
    if tools["prep_r"] and tools["pysh"]:
        try:
            run([tools["pysh"], tools["prep_r"], "-r", protein_pdb,
                 "-o", out_pdbqt, "-A", "checkhydrogens", "-U", "nphs_lps_waters"])
            if os.path.exists(out_pdbqt):
                if sanitize_pdbqt(out_pdbqt):        # fix long-name column overflow
                    log("صلّحت أعمدة receptor.pdbqt (اسم ذرة طويل)")
                ok = pdbqt_coords_ok(out_pdbqt)
        except Exception as e:
            log(f"prepare_receptor4 فشل ({e}) — obabel fallback")
    if not ok:                                       # robust fallback: clean columns
        run(f'obabel "{protein_pdb}" -O "{out_pdbqt}" -xr -p 7.4')
        sanitize_pdbqt(out_pdbqt)
    if not pdbqt_coords_ok(out_pdbqt):
        die("فشل تحضير الريسبتر (أعمدة الإحداثيات غير صالحة)")


# --------------------------- ligand prep ---------------------------------- #
def prepare_ligand(lig_pdb, out_pdbqt, tools, min_heavy=None):
    """PDB -> PDBQT preserving input coordinates, guaranteeing a single ROOT.
       Route: obabel->mol2 (perceive bonds + protonate) -> prepare_ligand4 ->
       validate (1 ROOT + enough atoms). Fallbacks: obabel direct, then rigid."""
    wd       = os.path.dirname(os.path.abspath(out_pdbqt)) or "."
    out_base = os.path.basename(out_pdbqt)
    mol2     = out_pdbqt.replace(".pdbqt", "_prep.mol2")

    # 1) clean single-molecule mol2 with bonds + hydrogens (uses CONECT from the
    #    input pdb, so complex crystal ligands don't fragment)
    mol2_ok = False
    if tools["obabel"]:
        try:
            run(f'obabel "{lig_pdb}" -O "{mol2}" -p 7.4 --partialcharge gasteiger', quiet=True)
            mol2_ok = os.path.exists(mol2)
        except Exception:
            mol2_ok = False
    src = os.path.basename(mol2) if mol2_ok else os.path.basename(lig_pdb)

    # 2) prepare_ligand4: flexible first, then rigid (-Z) which forces one ROOT
    if tools["prep_l"] and tools["pysh"]:
        for extra in ([], ["-Z"]):
            try:
                run([tools["pysh"], tools["prep_l"], "-l", src, "-o", out_base,
                     "-A", "hydrogens", "-U", "nphs"] + extra, cwd=wd)
            except Exception as e:
                log(f"prepare_ligand4 فشل ({e})")
            if pdbqt_ligand_ok(out_pdbqt, min_heavy):
                if extra: log("الليكاند حُضّر rigid (تعذّرت الشجرة المرنة)")
                return
    # 3) obabel direct fallback (from the CONECT-correct pdb -> single molecule)
    if tools["obabel"]:
        try:
            run(f'obabel "{lig_pdb}" -O "{out_pdbqt}" -p 7.4 --partialcharge gasteiger', quiet=True)
        except Exception:
            pass
        if pdbqt_ligand_ok(out_pdbqt, min_heavy): return
    die("فشل تحضير الليكاند (PDBQT غير صالح أو ناقص ذرات)")

def smiles_to_pdb(smiles, out_pdb):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.MolFromSmiles(smiles)
    if m is None: die(f"SMILES غلط: {smiles}")
    m = Chem.AddHs(m)
    if AllChem.EmbedMolecule(m, randomSeed=SEED) != 0:
        AllChem.EmbedMolecule(m, randomSeed=SEED, useRandomCoords=True)
    AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
    Chem.MolToPDBFile(m, out_pdb)


# --------------------------- docking -------------------------------------- #
def run_vina(receptor, ligand, center, size, out_pdbqt, log_txt, tools):
    cx, cy, cz = center; sx, sy, sz = size
    cmd = [tools["vina"],
           "--receptor", receptor, "--ligand", ligand,
           "--center_x", f"{cx:.3f}", "--center_y", f"{cy:.3f}", "--center_z", f"{cz:.3f}",
           "--size_x", f"{sx:.3f}", "--size_y", f"{sy:.3f}", "--size_z", f"{sz:.3f}",
           "--exhaustiveness", str(EXHAUSTIVENESS), "--num_modes", str(NUM_MODES),
           "--seed", str(SEED), "--out", out_pdbqt]
    out = run(cmd)
    with open(log_txt, "w") as f: f.write(out)
    # parse best affinity (mode 1)
    best = None
    for line in out.splitlines():
        p = line.split()
        if len(p) >= 2 and p[0] == "1":
            try: best = float(p[1]); break
            except: pass
    return best


def split_first_pose(dock_pdbqt, out_pdbqt):
    """Extract MODEL 1 from a multi-model vina output pdbqt."""
    lines, keep, started = [], [], False
    with open(dock_pdbqt) as f:
        for ln in f:
            if ln.startswith("MODEL"):
                if started: break
                started = True; continue
            if ln.startswith("ENDMDL"): break
            if started: keep.append(ln)
    with open(out_pdbqt, "w") as f: f.writelines(keep if keep else open(dock_pdbqt).readlines())


# --------------------------- RMSD (redock validation) --------------------- #
def _heavy_coords_pdbqt(path):
    """Heavy-atom coords in file order (AutoDock type in last column)."""
    xs = []
    for ln in open(path):
        if ln.startswith(("ATOM", "HETATM")):
            toks = ln.rstrip().split()
            atype = toks[-1] if toks else ""
            if atype.upper().startswith("H"):
                continue
            xs.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    return np.array(xs)

def compute_rmsd(native_pdbqt, pose_pdbqt, tools):
    """Redocking RMSD between the prepared native pose (Vina input) and the best
    docked pose (Vina output). Same atom order => direct heavy-atom RMSD (always
    finite). obrms is used only as a symmetry refinement if it returns smaller."""
    a = _heavy_coords_pdbqt(native_pdbqt)
    b = _heavy_coords_pdbqt(pose_pdbqt)
    n = min(len(a), len(b))
    if n == 0:
        return float("inf")
    order_rmsd = float(np.sqrt(((a[:n] - b[:n])**2).sum(1).mean()))
    if tools["obrms"] and tools["obabel"]:                # optional symmetry check
        try:
            ap = native_pdbqt.replace(".pdbqt", "_rref.pdb")
            bp = pose_pdbqt.replace(".pdbqt", "_rpose.pdb")
            run(f'obabel "{native_pdbqt}" -O "{ap}" -d', quiet=True)
            run(f'obabel "{pose_pdbqt}"  -O "{bp}" -d', quiet=True)
            v = float(run(f'obrms "{ap}" "{bp}"', quiet=True).strip().split()[-1])
            if v == v and v != float("inf") and v <= order_rmsd + 1e-6:
                return v
        except Exception:
            pass
    return order_rmsd


# --------------------------- rescoring & ligand efficiency ---------------- #
import re

def parse_score(out):
    """Extract the affinity (kcal/mol) from a vina --score_only run."""
    for key in ("Estimated Free Energy of Binding", "Affinity"):
        for ln in out.splitlines():
            if key in ln:
                m = re.search(r'[-+]?\d+\.\d+', ln)
                if m: return float(m.group())
    return None

def rescore(receptor, ligand_pdbqt, center, size, tools, scoring="vina"):
    """Re-score an existing pose in place (no docking) with vina or vinardo."""
    cx, cy, cz = center; sx, sy, sz = size
    cmd = [tools["vina"], "--receptor", receptor, "--ligand", ligand_pdbqt,
           "--center_x", f"{cx:.3f}", "--center_y", f"{cy:.3f}", "--center_z", f"{cz:.3f}",
           "--size_x", f"{sx:.3f}", "--size_y", f"{sy:.3f}", "--size_z", f"{sz:.3f}",
           "--scoring", scoring, "--score_only"]
    try:
        return parse_score(run(cmd, quiet=True))
    except Exception as e:
        log(f"rescore ({scoring}) فشل: {e}")
        return None

def ligand_efficiency(affinity, n_heavy):
    """LE = -ΔG / heavy-atom-count  (kcal/mol per heavy atom; higher = better)."""
    if affinity is None or not n_heavy: return None
    return -affinity / n_heavy

def n_heavy_from_smiles(smiles):
    from rdkit import Chem
    m = Chem.MolFromSmiles(smiles)
    return m.GetNumHeavyAtoms() if m else None

def fmt(x, nd=3):
    return f"{x:.{nd}f}" if isinstance(x, (int, float)) else "N/A"


# --------------------------- interaction analysis (PLIP) ------------------ #
def build_complex(receptor_pdb, ligand_pose_pdb, out_pdb):
    """Merge protein + docked ligand (as HETATM resname LIG) into one PDB."""
    prot = [l for l in open(receptor_pdb) if l.startswith(("ATOM","TER"))]
    lig  = []
    for l in open(ligand_pose_pdb):
        if l.startswith(("ATOM","HETATM")):
            l = "HETATM" + l[6:17] + "LIG Z 999" + l[26:]  # force resname/chain
            lig.append(l)
    with open(out_pdb, "w") as f:
        f.writelines(prot); f.write("TER\n"); f.writelines(lig); f.write("END\n")

def run_plip(complex_pdb, out_dir, tools):
    if not tools["plip"]:
        log("PLIP مو منصّب — تخطّي تحليل الارتباط"); return None
    os.makedirs(out_dir, exist_ok=True)
    try:
        run(f'plip -f "{complex_pdb}" -t -x -o "{out_dir}"')
    except Exception as e:
        log(f"PLIP فشل: {e}"); return None
    rep = os.path.join(out_dir, "report.txt")
    if os.path.exists(rep):
        log("---- بطاقات الارتباط (PLIP) ----")
        print(open(rep).read())
    return out_dir


# ============================== MAIN ====================================== #
def main():
    tools = discover_tools()
    pdb_id = PDB_ID.upper()
    outdir = os.path.join("results", pdb_id)
    os.makedirs(outdir, exist_ok=True)
    log(f"=== البروتين {pdb_id} → {outdir} ===")

    raw_pdb = os.path.join(outdir, f"{pdb_id}.pdb")
    download_pdb(pdb_id, raw_pdb)

    # 1) detect co-crystal ligand
    key, cands = detect_ligand(raw_pdb)
    log(f"مرشحين الكوكرستال: {cands}")
    log(f"الكوكرستال المختار: chain={key[0]} resname={key[1]} seq={key[2]}")

    # 2) split protein + native ligand
    protein_pdb = os.path.join(outdir, "protein.pdb")
    native_pdb  = os.path.join(outdir, "native_ligand.pdb")
    write_selection(raw_pdb, protein_pdb, ProteinSelect())
    extract_ligand_pdb(raw_pdb, key, native_pdb)

    # 3) receptor pdbqt
    receptor_pdbqt = os.path.join(outdir, "receptor.pdbqt")
    prepare_receptor(protein_pdb, receptor_pdbqt, tools)

    # 4) grid box from native ligand
    center, size = grid_box_from_ligand(native_pdb, BOX_PADDING)
    with open(os.path.join(outdir, "grid_box.txt"), "w") as f:
        f.write(f"center_x {center[0]:.3f}\ncenter_y {center[1]:.3f}\ncenter_z {center[2]:.3f}\n"
                f"size_x {size[0]:.3f}\nsize_y {size[1]:.3f}\nsize_z {size[2]:.3f}\n")
    log(f"BOX center={np.round(center,2).tolist()} size={np.round(size,2).tolist()}")

    # 5) REDOCK native ligand → RMSD validation
    log("=== REDOCK validation ===")
    native_pdbqt = os.path.join(outdir, "native_ligand.pdbqt")
    prepare_ligand(native_pdb, native_pdbqt, tools, min_heavy=cands[key])
    redock_out = os.path.join(outdir, "redock_out.pdbqt")
    redock_aff = run_vina(receptor_pdbqt, native_pdbqt, center, size, redock_out,
                          os.path.join(outdir, "redock_log.txt"), tools)
    pose1_pdbqt = os.path.join(outdir, "redock_pose1.pdbqt")
    pose1_pdb   = os.path.join(outdir, "redock_pose1.pdb")
    split_first_pose(redock_out, pose1_pdbqt)
    run(f'obabel "{pose1_pdbqt}" -O "{pose1_pdb}"', quiet=True)
    rmsd = compute_rmsd(native_pdbqt, pose1_pdbqt, tools)
    verdict = "موثوق ✓" if rmsd <= RMSD_CUTOFF else "ضعيف ✗ (راجع البوكس/التحضير)"
    log(f"REDOCK RMSD = {rmsd:.2f} Å  →  {verdict}")
    # rescoring + ligand efficiency for the CO-CRYSTAL (for fair comparison)
    nH_native   = cands[key]
    nat_vina    = rescore(receptor_pdbqt, pose1_pdbqt, center, size, tools, "vina")
    nat_vinardo = rescore(receptor_pdbqt, pose1_pdbqt, center, size, tools, "vinardo")
    le_native   = ligand_efficiency(redock_aff, nH_native)
    log(f"الكوكرستال {key[1]}: dock={redock_aff}  vina_rescore={nat_vina}  "
        f"vinardo={nat_vinardo}  LE={fmt(le_native)}")

    # 6) DOCK query compound
    log("=== docking المركب ===")
    query_pdb   = os.path.join(outdir, "query.pdb")
    query_pdbqt = os.path.join(outdir, "query.pdbqt")
    smiles_to_pdb(QUERY_SMILES, query_pdb)
    prepare_ligand(query_pdb, query_pdbqt, tools, min_heavy=n_heavy_from_smiles(QUERY_SMILES))
    query_out = os.path.join(outdir, "query_docked.pdbqt")
    best_aff = run_vina(receptor_pdbqt, query_pdbqt, center, size, query_out,
                        os.path.join(outdir, "query_log.txt"), tools)
    log(f"افضل affinity للمركب = {best_aff} kcal/mol")

    # 7) interaction analysis on best pose
    qpose_pdbqt = os.path.join(outdir, "query_pose1.pdbqt")
    qpose_pdb   = os.path.join(outdir, "query_pose1.pdb")
    split_first_pose(query_out, qpose_pdbqt)
    run(f'obabel "{qpose_pdbqt}" -O "{qpose_pdb}"', quiet=True)
    # rescoring + ligand efficiency for the QUERY compound
    nH_query   = n_heavy_from_smiles(QUERY_SMILES)
    q_vina     = rescore(receptor_pdbqt, qpose_pdbqt, center, size, tools, "vina")
    q_vinardo  = rescore(receptor_pdbqt, qpose_pdbqt, center, size, tools, "vinardo")
    le_query   = ligand_efficiency(best_aff, nH_query)
    log(f"المركب: dock={best_aff}  vina_rescore={q_vina}  "
        f"vinardo={q_vinardo}  LE={fmt(le_query)}")
    complex_pdb = os.path.join(outdir, "complex.pdb")
    build_complex(protein_pdb, qpose_pdb, complex_pdb)
    run_plip(complex_pdb, os.path.join(outdir, "interactions"), tools)

    # 8) summary + comparison table + CSV
    summary = os.path.join(outdir, "SUMMARY.txt")
    q_name  = "QUERY"
    c_name  = f"CO-CRYSTAL({key[1]})"
    rows = [
        ("Heavy atoms",                 str(nH_query),  str(nH_native)),
        ("Vina affinity (docking)",     fmt(best_aff),  fmt(redock_aff)),
        ("Vina rescore (score_only)",   fmt(q_vina),    fmt(nat_vina)),
        ("Vinardo rescore",             fmt(q_vinardo), fmt(nat_vinardo)),
        ("Ligand Efficiency (-dG/HAC)", fmt(le_query),  fmt(le_native)),
        ("Redock RMSD (A)",             "-",            fmt(rmsd, 2)),
    ]
    tbl = f"{'':<30}{q_name:>18}{c_name:>22}\n"
    for lbl, a, b in rows:
        tbl += f"{lbl:<30}{a:>18}{b:>22}\n"
    with open(summary, "w") as f:
        f.write(textwrap.dedent(f"""\
        ================ DOCKING SUMMARY : {pdb_id} ================
        Co-crystal ligand : chain {key[0]}  resname {key[1]}  seq {key[2]}
        Grid center (x,y,z): {center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}
        Grid size   (x,y,z): {size[0]:.3f}, {size[1]:.3f}, {size[2]:.3f}
        Redock validation : {rmsd:.2f} A   ({verdict})
        Query SMILES      : {QUERY_SMILES}
        ========================== COMPARISON ==========================
        """))
        f.write(tbl)
        f.write(textwrap.dedent("""\
        ================================================================
        Notes:
         - Vina affinity  = best docked mode (mode 1) energy, kept as-is.
         - Rescore        = same pose re-evaluated in place (--score_only).
         - Ligand Eff.    = -affinity / heavy-atom-count (higher = better);
                            fairer than raw affinity when sizes differ.
        Files: receptor.pdbqt, query_docked.pdbqt, complex.pdb, scores.csv,
               interactions/report.txt (بطاقات الارتباط)
        ================================================================
        """))
    # CSV (accumulate rows across proteins by hand, or append externally)
    csv_path = os.path.join(outdir, "scores.csv")
    with open(csv_path, "w") as f:
        f.write("pdb_id,compound,resname,heavy_atoms,vina_dock,vina_rescore,"
                "vinardo_rescore,ligand_efficiency,redock_rmsd\n")
        f.write(f"{pdb_id},query,-,{nH_query},{fmt(best_aff)},{fmt(q_vina)},"
                f"{fmt(q_vinardo)},{fmt(le_query)},\n")
        f.write(f"{pdb_id},cocrystal,{key[1]},{nH_native},{fmt(redock_aff)},"
                f"{fmt(nat_vina)},{fmt(nat_vinardo)},{fmt(le_native)},{fmt(rmsd,2)}\n")
    log("=== خلصت ✓ ===  شوف: " + summary + "  و  " + csv_path)
    print(open(summary).read())


if __name__ == "__main__":
    main()
