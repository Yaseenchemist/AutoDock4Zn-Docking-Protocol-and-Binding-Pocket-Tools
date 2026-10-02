# بروتوكول AutoDock4Zn — النسخة النهائية (تركيب محفوظ)
### إرساء إنزيمات الزنك — تحضير بـ Open Babel + استرجاع الآصرة الصحيحة للـ poses

> **جديد ومهم:** كل الليكندات (كوكرستال + جديدة) تُحضّر بـ **Open Babel فقط** (لا `prepare_ligand4`).
> بعد الدوكنك، كل pose يُعاد بناؤه بالآصرة الصحيحة من ملف الأصل (template) → **لا يبقى أي ملف بتركيب متلف**.
> السبب: PDBQT/PDB لا يخزّنان الآصرة، فالعارض (MOE) يخمّنها غلط ويفتح الحلقات. الحل: نرجّع الآصرة من الأصل.

---

# ★ القسم 0 — صندوق الإعدادات

```bash
cd ~/Desktop/shereen          # ← (1) مجلد العمل
PROTEIN=protein.pdb           # ← (2) البروتين الخام (زنك + كوكرستال)
COCRYSTAL=SHH                 # ← (3) كود الكوكرستال (3 أحرف)
LIGSDF=14lig.sdf              # ← (4) ملف الليكندات SDF
NPTS=50                       # ← (5) حجم الصندوق (50 = ~19Å)
NPOSES_VAL=50                 # ← (6) poses الريدوك/التحقق (يُنصح يبقى 50 — دقة)
NPOSES_DOCK=10                # ← (7) poses المركبات الجديدة (10 = أسرع)

conda activate ad4zn
cp "$AD4ZN"/AD4Zn.dat .

# بايثون فيه RDKit (لاسترجاع الآصرة) — عدّل المسار لمسارك
PYR=/home/yaseen/SOFTWARE/miniforge3/bin/python3
```

---

# ⚙️ القسم 0.A — تعديل عتبة الزنك (اختياري، قبل الريدوك)

> افتراضياً لا تعدّل (2.5 Å صحيح). عدّل فقط لو تتوقع دونر على 2.5–3.0 Å. الوقت: ثوانٍ + إعادة القسم 1.

```bash
sed -i 's/cutoff = 2.5/cutoff = 2.8/' ~/SOFTWARE/AD4Zn/zinc_pseudo.py   # (1) الأهم — يغيّر الدوكنك
sed -i 's/if d<3.6/if d<5.0/'        /tmp/zn_shell.py                   # (2) عرض فقط
sed -i 's/d<=2.5/d<=2.8/g'           /tmp/zn_shell.py                   # (3) فحص فقط
```

---

# القسم 0.B — السكربتات المساعدة (تُكتب مرّة)

**فحص الزنك:**
```bash
cat > /tmp/zn_shell.py << 'EOF'
import sys, math
f = sys.argv[1]; atoms=[]; zn=None
for L in open(f):
    if L[:6] in ("ATOM  ","HETATM"):
        try: x=float(L[30:38]); y=float(L[38:46]); z=float(L[46:54])
        except ValueError: continue
        name=L[12:16].strip(); res=L[17:20].strip(); seq=L[22:26].strip(); t=L.split()[-1]
        rec=[name,res,seq,t,x,y,z]; atoms.append(rec)
        if t.upper()=="ZN" or name.upper()=="ZN": zn=rec
if not zn: print("No ZN found in",f); sys.exit()
zx,zy,zz=zn[4],zn[5],zn[6]; near=[]
for a in atoms:
    if a is zn: continue
    d=math.sqrt((a[4]-zx)**2+(a[5]-zy)**2+(a[6]-zz)**2)
    if d<3.6: near.append((d,a))
near.sort(); valid={'O','OA','NA','N','S','SA','MG','MN','CA','FE','CU'}
print("Zn at (%.3f %.3f %.3f)\n"%(zx,zy,zz)); print("dist   atom  res   seq  type"); n25=0
for d,a in near:
    don=a[3] in valid
    flag="  <-- DONOR (<=2.5)" if(don and d<=2.5) else("  (donor, 2.5-3.6)" if don else "")
    print("%5.2f  %-4s  %-4s  %-4s %-3s%s"%(d,a[0],a[1],a[2],a[3],flag))
    if don and d<=2.5: n25+=1
print("\nValid donors within 2.5 A = %d"%n25)
EOF
```

**استرجاع الآصرة لأي pose (template + pose → SDF/MOL2/PDB صحيح):**
```bash
cat > ~/SOFTWARE/AD4Zn/rebuild_pose.py << 'EOF'
#!/usr/bin/env python3
import sys, subprocess
from rdkit import Chem
from rdkit.Geometry import Point3D
template_sdf, pose_pdbqt, out_prefix = sys.argv[1], sys.argv[2], sys.argv[3]
subprocess.run(["obabel", pose_pdbqt, "-O", "/tmp/_pose_tmp.pdb"], check=True,
               stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
tmpl = Chem.RemoveHs(Chem.MolFromMolFile(template_sdf))
pose = Chem.MolFromPDBFile("/tmp/_pose_tmp.pdb", removeHs=True, sanitize=False)
if tmpl is None or pose is None: sys.exit("ERROR reading files")
if tmpl.GetNumAtoms() != pose.GetNumAtoms():
    sys.exit("ATOM MISMATCH %d vs %d"%(tmpl.GetNumAtoms(), pose.GetNumAtoms()))
conf = Chem.Conformer(tmpl.GetNumAtoms()); pc = pose.GetConformer()
for i in range(tmpl.GetNumAtoms()):
    p = pc.GetAtomPosition(i); conf.SetAtomPosition(i, Point3D(p.x, p.y, p.z))
tmpl.RemoveAllConformers(); tmpl.AddConformer(conf, assignId=True)
Chem.SanitizeMol(tmpl)
Chem.MolToMolFile(tmpl, out_prefix+".sdf")
for fmt in ("mol2","pdb"):
    subprocess.run(["obabel", out_prefix+".sdf", "-O", out_prefix+"."+fmt], check=True,
                   stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
print("OK", out_prefix, "rings:", tmpl.GetRingInfo().NumRings())
EOF
```

**دمج الكومبلكس مع CONECT (بروتين+زنك+ليكند صحيح):**
```bash
cat > ~/SOFTWARE/AD4Zn/merge_complex.py << 'EOF'
#!/usr/bin/env python3
import sys
receptor_pdb, ligand_pdb, out_pdb = sys.argv[1], sys.argv[2], sys.argv[3]
rec=[]; maxser=0
for L in open(receptor_pdb):
    if L[:6] in ("ATOM  ","HETATM"):
        rec.append(L.rstrip("\n"))
        try: maxser=max(maxser,int(L[6:11]))
        except: pass
lig=[]; con=[]
for L in open(ligand_pdb):
    if L[:6] in ("ATOM  ","HETATM"):
        s=int(L[6:11])+maxser
        lig.append("HETATM%5d%s"%(s, L[11:].rstrip("\n")))
    elif L[:6]=="CONECT":
        nums=[int(L[i:i+5])+maxser for i in range(6,len(L.rstrip()),5) if L[i:i+5].strip()]
        con.append("CONECT"+"".join("%5d"%n for n in nums))
with open(out_pdb,"w") as o:
    for L in rec: o.write(L+"\n")
    o.write("TER\n")
    for L in lig: o.write(L+"\n")
    for L in con: o.write(L+"\n")
    o.write("END\n")
print("OK", out_pdb, "| ligand atoms:", len(lig), "| CONECT:", len(con))
EOF
```

---

## القسم 1 — تحضير البروتين + TZ

```bash
grep "$COCRYSTAL" $PROTEIN > cocrystal.pdb
grep -v "$COCRYSTAL" $PROTEIN > protein_clean.pdb
prep_rec -r protein_clean.pdb -o protein.pdbqt -A hydrogens
python3 "$AD4ZN"/zinc_pseudo.py -r protein.pdbqt -o protein_tz.pdbqt
grep -E "ZN|TZ" protein_tz.pdbqt
python3 /tmp/zn_shell.py protein.pdbqt          # لازم 3 دونرز
```

---

## القسم 2 — تحديد البوكت (تلقائي على الزنك)

```bash
read CX CY CZ <<< $(grep -E "ZN .* ZN " protein_tz.pdbqt | awk '{print $7, $8, $9}')
echo "gridcenter = $CX $CY $CZ | box = ${NPTS}pts"
```

---

## القسم 3 — تحضير الليكندات بـ Open Babel + حفظ الأصل (template)

> **كله Open Babel** — لا `prepare_ligand4`. نحفظ الأصل لاسترجاع الآصرة لاحقاً.

```bash
# pdbqt للدوكنك (Open Babel)
obabel $LIGSDF -O lig_.pdbqt -m --partialcharge gasteiger -p 7.4

# الأصل كـ template (نفس ترتيب الذرات — مهم للاسترجاع)
obabel $LIGSDF -O lig_orig_.sdf -m

NLIG=$(ls lig_*.pdbqt | wc -l); echo "عدد الليكندات = $NLIG"
LT=$(grep -hE "^ATOM|^HETATM" lig_*.pdbqt | awk '{print $NF}' | sort -u | paste -sd, -)
echo "ligand types = [$LT]"

mkdir -p Ligands_pdbqt && cp lig_*.pdbqt Ligands_pdbqt/
```

---

## القسم 4 — GPF + maps

```bash
prep_gpf_zn -l lig_1.pdbqt -r protein_tz.pdbqt \
            -p gridcenter=$CX,$CY,$CZ \
            -p npts=$NPTS,$NPTS,$NPTS -p spacing=0.375 \
            -p ligand_types="$LT" -o protein_tz.gpf
autogrid4 -p protein_tz.gpf -l protein_tz.glg
tail -2 protein_tz.glg
```

---

## القسم 5 — التحقق (Redock) — كوكرستال بـ Open Babel

```bash
# تحضير الكوكرستال بـ Open Babel (لا prepare_ligand4) + template
obabel cocrystal.pdb -O cocrystal.pdbqt --partialcharge gasteiger
obabel cocrystal.pdb -O cocrystal_orig.sdf

prep_dpf -l cocrystal.pdbqt -r protein_tz.pdbqt -o redock.dpf
sed -i "s/^ga_run.*/ga_run $NPOSES_VAL/" redock.dpf
autodock4 -p redock.dpf -l redock.dlg

# جدول: Pose | Energy | RMSD | Zn-dist
awk '/^DOCKED/{sub(/^DOCKED: /,""); print}' redock.dlg > all_docked.pdbqt
csplit -sz -f _pose_ -b "%03d.pdbqt" all_docked.pdbqt '/^MODEL/' '{*}' 2>/dev/null
grep -E "ZN .* ZN " protein_tz.pdbqt > zn_only.pdb
echo "Pose | Energy | RMSD(Å) | Zn-dist(Å)"; n=0
for f in _pose_*.pdbqt; do
  n=$((n+1))
  E=$(grep "Estimated Free Energy" "$f" | awk '{print $8}' | head -1)
  R=$(obrms cocrystal.pdbqt "$f" 2>/dev/null | awk '{print $NF}')
  cat zn_only.pdb "$f" > /tmp/zp.pdbqt
  Z=$(python3 /tmp/zn_shell.py /tmp/zp.pdbqt 2>/dev/null | awk '/^ *[0-9]+\.[0-9]+ /{print $1; exit}')
  printf "%-4s | %-7s | %-7s | %s\n" "$n" "$E" "$R" "$Z"
done
rm -f _pose_*.pdbqt
```
> التحقق ناجح: أي pose RMSD < 2 Å أو يحافظ على تنسيق الزنك (~2.0–2.5 Å).

---

## القسم 6 — دوكنك الـ 14

```bash
for i in $(seq 1 $NLIG); do
  echo ">>> [$( date +%H:%M )] lig_$i ..."
  prep_dpf -l lig_$i.pdbqt -r protein_tz.pdbqt -o dock_lig_$i.dpf
  sed -i "s/^ga_run.*/ga_run $NPOSES_DOCK/" dock_lig_$i.dpf
  autodock4 -p dock_lig_$i.dpf -l dock_lig_$i.dlg
done
echo "===== خلصوا ====="
```

---

## القسم 7 — جدول الطاقات + الزنك (مرتّب بالـ affinity)

```bash
cat > /tmp/scan_all.py << EOF
import math, glob, re
ZN=($CX,$CY,$CZ); DONOR={'N','NA','OA','O','SA','S'}
def parse(path):
    poses=[]; atoms=[]; E=None
    for line in open(path):
        if line[:7]!="DOCKED:": continue
        c=line[8:] if line[:8]=="DOCKED: " else line[7:]
        if c.startswith("MODEL"): atoms=[]; E=None
        elif "Estimated Free Energy of Binding" in c:
            m=re.search(r'=\s*([-0-9.]+)\s*kcal',c);  E=float(m.group(1)) if m else E
        elif c[:4]=="ATOM" or c[:6]=="HETATM":
            try:
                x=float(c[30:38]);y=float(c[38:46]);z=float(c[46:54]);t=c.split()[-1]
                atoms.append((t,x,y,z))
            except: pass
        elif c.startswith("ENDMDL"):
            best=None
            for t,x,y,z in atoms:
                if t in DONOR:
                    d=math.sqrt((x-ZN[0])**2+(y-ZN[1])**2+(z-ZN[2])**2)
                    if best is None or d<best: best=d
            poses.append((E,best))
    return poses
rows=[]
for p in glob.glob("dock_lig_*.dlg"):
    ps=parse(p); n=int(re.search(r'(\d+)',p).group(1))
    valid=[x for x in ps if x[0] is not None]
    if not valid: continue
    r1=min(valid,key=lambda x:x[0]); coord=[x for x in ps if x[1] is not None]
    bc=min(coord,key=lambda x:x[1]) if coord else (None,None)
    rows.append((n,r1[0],r1[1],bc[1]))
rows.sort(key=lambda r:r[1])
print("Rank | Ligand | BestE | rank1_Zn | BEST_Zn | coord?")
for i,(n,E,z1,zb) in enumerate(rows,1):
    yn="YES" if (zb and zb<=2.5) else "no"
    print("%-4d | lig_%-3d | %-7.2f | %-8s | %-7s | %s"%(i,n,E,
        "%.2f"%z1 if z1 else "-","%.2f"%zb if zb else "-",yn))
EOF
python3 /tmp/scan_all.py | tee results_summary.txt
```

---

## القسم 8 — بناء `AD4Zn_Results` (تركيب صحيح + complex مع CONECT)

> لكل ليكند: نستخرج أوطأ pose (MODEL 1) → نعيد بناءه بالآصرة الصحيحة من `lig_orig_N.sdf` → نبني complex مع CONECT.

```bash
cd ~/Desktop/shereen

# ========== 1) استبدال rebuild_pose.py بالنسخة المصححة ==========
cat > ~/SOFTWARE/AD4Zn/rebuild_pose.py << 'EOF'
#!/usr/bin/env python3
import sys, subprocess
from rdkit import Chem
from rdkit.Geometry import Point3D
template_sdf, pose_pdbqt, out_prefix = sys.argv[1], sys.argv[2], sys.argv[3]
subprocess.run(["obabel", pose_pdbqt, "-O", "/tmp/_pose_noH.pdb", "-d"],
               check=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
pose = Chem.MolFromPDBFile("/tmp/_pose_noH.pdb", removeHs=True, sanitize=False)
tmpl = Chem.RemoveHs(Chem.MolFromMolFile(template_sdf))
if tmpl is None or pose is None: sys.exit("ERROR reading template/pose")
t_heavy=[a.GetIdx() for a in tmpl.GetAtoms() if a.GetAtomicNum()>1]
p_heavy=[a.GetIdx() for a in pose.GetAtoms() if a.GetAtomicNum()>1]
if len(t_heavy)!=len(p_heavy):
    sys.exit("ATOM MISMATCH heavy %d vs %d"%(len(t_heavy),len(p_heavy)))
conf=Chem.Conformer(tmpl.GetNumAtoms()); pc=pose.GetConformer()
for ti,pi in zip(t_heavy,p_heavy):
    p=pc.GetAtomPosition(pi); conf.SetAtomPosition(ti,Point3D(p.x,p.y,p.z))
tmpl.RemoveAllConformers(); tmpl.AddConformer(conf,assignId=True)
Chem.SanitizeMol(tmpl)
Chem.MolToMolFile(tmpl, out_prefix+".sdf")
for fmt in ("mol2","pdb"):
    subprocess.run(["obabel", out_prefix+".sdf", "-O", out_prefix+"."+fmt],
                   check=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
print("OK", out_prefix, "| heavy:", len(t_heavy), "| rings:", tmpl.GetRingInfo().NumRings())
EOF

# ========== 2) إعادة بناء AD4Zn_Results ==========
mkdir -p AD4Zn_Results/_Receptor
cp protein_tz.pdbqt AD4Zn_Results/_Receptor/
grep -vE "TZ" protein_tz.pdbqt | cut -c1-66 | grep -E "^ATOM|^HETATM" > AD4Zn_Results/_Receptor/receptor_Zn.pdb
cp AD4Zn.dat protein_tz.gpf protein_tz.maps.fld AD4Zn_Results/_Receptor/ 2>/dev/null
grep -E "ZN .* ZN " protein_tz.pdbqt > zn_only.pdb
REC=AD4Zn_Results/_Receptor/receptor_Zn.pdb

for i in $(seq 1 $NLIG); do
  dlg=dock_lig_$i.dlg; [ -f "$dlg" ] || continue
  D=AD4Zn_Results/lig_$i; mkdir -p "$D"
  cp lig_$i.pdbqt lig_orig_$i.sdf dock_lig_$i.dpf dock_lig_$i.dlg "$D/"
  awk '/^DOCKED/{sub(/^DOCKED: /,""); print}' "$dlg" > "$D/all_poses.pdbqt"
  csplit -sz -f "$D/_p_" -b "%03d.pdbqt" "$D/all_poses.pdbqt" '/^MODEL/' '{*}' 2>/dev/null
  mv "$D/_p_000.pdbqt" "$D/best_pose_raw.pdbqt"; rm -f "$D/_p_"*.pdbqt
  $PYR ~/SOFTWARE/AD4Zn/rebuild_pose.py lig_orig_$i.sdf "$D/best_pose_raw.pdbqt" "$D/best_pose"
  $PYR ~/SOFTWARE/AD4Zn/merge_complex.py $REC "$D/best_pose.pdb" "$D/complex_full.pdb"
  rm -f "$D/best_pose_raw.pdbqt"
  E=$(grep "Estimated Free Energy" "$dlg" | awk '{print $8}' | sort -n | head -1)
  cat zn_only.pdb "$D/best_pose.pdb" > /tmp/zp.pdbqt
  Dz=$($PYR /tmp/zn_shell.py /tmp/zp.pdbqt 2>/dev/null | awk '/^ *[0-9]+\.[0-9]+ /{print $1; exit}')
  printf "lig_%s\nBest Energy (affinity) = %s kcal/mol\nZn-closest = %s A\n" "$i" "$E" "$Dz" > "$D/INFO.txt"
done
echo "===== تم. محتوى lig_1: ====="; ls -1 AD4Zn_Results/lig_1/

> كل فولدر فيه: `best_pose.sdf/.mol2/.pdb` (تركيب صحيح) + `complex_full.pdb` (مع CONECT) + dlg + INFO.

---

## القسم 9 — بناء `Cocrystal_Redock_Results` (تركيب صحيح)

```bash
mkdir -p Cocrystal_Redock_Results
awk '/^DOCKED/{sub(/^DOCKED: /,""); print}' redock.dlg > all_docked.pdbqt
csplit -sz -f _cp_ -b "%03d.pdbqt" all_docked.pdbqt '/^MODEL/' '{*}' 2>/dev/null
mv _cp_000.pdbqt best_redock_raw.pdbqt; rm -f _cp_*.pdbqt

# استرجاع الآصرة + complex مع CONECT
$PYR ~/SOFTWARE/AD4Zn/rebuild_pose.py cocrystal_orig.sdf best_redock_raw.pdbqt Cocrystal_Redock_Results/best_pose
$PYR ~/SOFTWARE/AD4Zn/merge_complex.py AD4Zn_Results/_Receptor/receptor_Zn.pdb \
     Cocrystal_Redock_Results/best_pose.pdb Cocrystal_Redock_Results/complex_full.pdb
cp cocrystal.pdb cocrystal_orig.sdf redock.dpf redock.dlg Cocrystal_Redock_Results/
rm -f best_redock_raw.pdbqt
ls -lh Cocrystal_Redock_Results/
```

> الكومبلكس = أوطأ طاقة. لو تريد الـ pose المتنسّق بدلاً منه، بدّل `_cp_000` برقم الـ pose المتنسّق من جدول القسم 5.

---

## القسم 9.5 — تنظيم الملفات بفولدرات (تقليل الضوضاء) 🗂️

> بدل ما تبقى الملفات مبعثرة بالجذر، ننقلها لفولدرات منظّمة. **النتائج النهائية بـ AD4Zn_Results / Cocrystal_Redock_Results تبقى بمكانها.**

```bash
cd ~/Desktop/shereen

# (1) فولدرات التنظيم
mkdir -p _Inputs _Receptor_prep _Ligands_prep _Raw_docking _Temp _Analysis

# (2) المدخلات الأصلية
mv "$PROTEIN" protein_clean.pdb cocrystal.pdb cocrystal_orig.sdf $LIGSDF _Inputs/ 2>/dev/null
mv *.smi smiles_*.txt _Inputs/ 2>/dev/null

# (3) تحضير البروتين (الزنك + maps + TZ)
mv protein.pdbqt protein_tz.pdbqt protein_tz.gpf protein_tz.glg \
   protein_tz.*.map protein_tz.maps.fld protein_tz.maps.xyz AD4Zn.dat _Receptor_prep/ 2>/dev/null

# (4) الليكندات المحضّرة + الأصول
mv lig_*.pdbqt lig_orig_*.sdf cocrystal.pdbqt Ligands_pdbqt _Ligands_prep/ 2>/dev/null

# (5) ملفات الدوكنك الخام (dpf + dlg)
mv dock_lig_*.dpf dock_lig_*.dlg redock.dpf redock.dlg _Raw_docking/ 2>/dev/null

# (6) ملفات التحليل
mv results_summary.txt zn_scan.txt _Analysis/ 2>/dev/null

# (7) المؤقت — يُحذف نهائياً (ما إله قيمة)
rm -f pose_*.pdbqt zp_*.pdbqt all_docked.pdbqt _pose_*.pdbqt _cp_*.pdbqt \
      *_raw.pdbqt zn_only.pdb protein_clean.pdb /tmp/zp.pdbqt /tmp/_pose_tmp.pdb
rmdir _Temp 2>/dev/null

echo "===== الجذر بعد التنظيم ====="
ls -1
```

**الهيكل النهائي للجذر (نظيف):**
```
~/Desktop/shereen/
├── AD4Zn_Results/            ← النتائج (فولدر لكل ليكند + بروتين)
├── Cocrystal_Redock_Results/ ← نتائج التحقق
├── _Inputs/                  ← المدخلات الأصلية
├── _Receptor_prep/           ← البروتين + الزنك + maps
├── _Ligands_prep/            ← الليكندات المحضّرة + الأصول
├── _Raw_docking/             ← dpf + dlg الخام
└── _Analysis/                ← الجداول والملخصات
```

> الفولدرات تبدأ بـ `_` حتى تتجمّع بالأعلى ومرتّبة. النتائج المهمة بدون `_` تبقى بارزة.

---

## القسم 10 — عرض الـ Interactions (Discovery Studio / MOE)

> الآن `complex_full.pdb` فيه CONECT → التركيب يظهر صحيح (حلقات سليمة).

1. افتح `AD4Zn_Results/lig_X/complex_full.pdb`
2. الليكند (HETATM) → **Define Ligand**
3. **Receptor-Ligand Interactions → 2D Diagram**

> بديل أنظف: افتح `_Receptor/receptor_Zn.pdb` + `lig_X/best_pose.sdf` (SDF فيه آصرة صريحة دائماً).

---

## ملاحظات

- **تركيب الـ pose المتلف** سببه فقدان الآصرة بالـ PDBQT/PDB. الحل (هذا البروتوكول): استرجاع الآصرة من الأصل + CONECT. الذرات كانت دائماً بأماكنها الصحيحة.
- **كل تحضير الليكندات Open Babel** (لا `prepare_ligand4` الذي يكسر الحلقات أحياناً).
- **التنسيق الضعيف** لمركبات thione/hydrazone متوقع كيميائياً (ليست hydroxamate).
- **التذبذب** بين تشغيلتين (RMSD) طبيعي (seed عشوائي).

## استكشاف أخطاء

| العَرَض | الحل |
|---|---|
| `ATOM MISMATCH` بالـ rebuild | الـ template (lig_orig) من مصدر مختلف عن الـ pdbqt — استعملهم من نفس SDF |
| تركيب متلف بـ MOE | استعمل `complex_full.pdb` الجديد (مع CONECT) أو `best_pose.sdf` |
| `Wrote 0 TZ atoms` | احذف الكوكرستال / دونرز ≠ 3 |
| RMSD فاضي لـ 10+ | csplit (موجود) مو awk |
| `autogrid4 not found` | `conda activate ad4zn` |

---
*بروتوكول نهائي — Open Babel فقط + استرجاع الآصرة + complex مع CONECT.*
