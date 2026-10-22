"""Rebuild the four main and fourteen supplementary figures, without new inference."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parents[1]
STAGES={
 'flow': "import dataset_flow; dataset_flow.build(); import legacy_figures as m; m.figure1()",
 'legacy': "import legacy_figures as m; old=m.save; m.save=lambda f,n,*a:old(f,{'figure3':'figureS11','figure4':'figureS14'}[n],*a); m.figure3(); m.figure4()",
 'supplement': "import supplement_figures",
 'genomovars': "import v5_figures as m; old=m.save; m.save=lambda f,n,*a:old(f,{'figureS9':'figureS10','figureS10':'figureS3'}[n],*a)\nwith m.plt.rc_context({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'}):\n a=m.genomovars(); b=m.subspecies()\n(m.RUN/'checks/genomovars.json').write_text(m.json.dumps({'genomovars':a,'subspecies':b},indent=2))",
 'trees': "import runpy; runpy.run_module('v7_treefigure',run_name='__main__')",
 'kraken': "import v6_figures as m; a=m.kraken(); old=m.save; m.save=lambda f,n:old(f,'figureS12'); m.diagram(); (m.RUN/'checks/kraken.json').write_text(m.json.dumps(a,indent=2))",
 'classification': "import v7_figure4 as m; a=m.pterobacter(); (m.RUN/'checks/classification.json').write_text(m.json.dumps(a,indent=2))",
}
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True,help='A new directory for generated files')
    args=parser.parse_args(); out=args.output.resolve()
    if out.exists():parser.error('Output must not exist; retained runs are never overwritten')
    out.mkdir(parents=True)
    for name in ['outputs/tables','outputs/figures/source-data','checks','logs','inputs']:(out/name).mkdir(parents=True)
    for folder,target in [('tables','outputs/tables'),('source-data','outputs/figures/source-data')]:
        for source in (PROJECT/'data/figure-seeds'/folder).iterdir():shutil.copyfile(source,out/target/source.name)
    env=os.environ.copy();env.update(ENTEROBACTER_REBUILD_OUTPUT=str(out),MPLCONFIGDIR=str(out/'outputs/mpl-cache'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
    for name,code in STAGES.items():
        print('Rebuilding',name,flush=True)
        with (out/'logs'/f'{name}.log').open('w') as log:
            subprocess.run([sys.executable,'-c',code],cwd=HERE,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
        if name=='supplement':
            mapping={1:1,2:4,3:5,4:9,5:8,6:7,7:6,8:2}
            figs=out/'outputs/figures'
            for old,new in mapping.items():
                for ext in ['png','svg','pdf']:(figs/f'figureS{old}.{ext}').rename(figs/f'pendingS{new}.{ext}')
            for new in mapping.values():
                for ext in ['png','svg','pdf']:(figs/f'pendingS{new}.{ext}').rename(figs/f'figureS{new}.{ext}')
    import numpy as np
    from PIL import Image
    checks=[]
    for stem in ['figure'+str(i) for i in range(1,5)]+['figureS'+str(i) for i in range(1,15)]:
        actual=out/'outputs/figures'/f'{stem}.png';expected=PROJECT/'expected/figures'/f'{stem}.png'
        a=np.asarray(Image.open(actual));b=np.asarray(Image.open(expected))
        checks.append({'figure':stem,'pixels_identical':a.shape==b.shape and bool(np.array_equal(a,b)), 'generated_sha256':hashlib.sha256(actual.read_bytes()).hexdigest(),'expected_sha256':hashlib.sha256(expected.read_bytes()).hexdigest()})
    result={'figures':checks,'all_pixels_identical':all(x['pixels_identical'] for x in checks),'python':sys.version}
    (out/'checks/figure-reproduction.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    if not result['all_pixels_identical']:raise SystemExit('Figure mismatch: inspect checks/figure-reproduction.json')
if __name__=='__main__':main()
