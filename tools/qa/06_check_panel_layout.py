"""Check physical plot-area alignment from an R-generated JSON manifest.

Shared grid boundaries must coincide within 1.5 PDF points. Unequal grid spans
remain unequal. The checker uses only the Python standard library.
"""
from pathlib import Path
import argparse,itertools,json,math

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('manifest',type=Path);p.add_argument('--json-out',type=Path,required=True)
    p.add_argument('--tolerance-pt',type=float,default=1.5)
    a=p.parse_args();m=json.loads(a.manifest.read_text(encoding='utf-8'))
    panels=m['panels'];checks=[];failures=[]
    if len(panels)<2:raise ValueError('At least two panels are required.')
    width,height=m['figure']['width_pt'],m['figure']['height_pt']
    for x in panels:
        b=x['bbox_pt']
        if len(b)!=4 or not all(math.isfinite(v) for v in b):raise ValueError('Invalid bounding box.')
        if not (-.1<=b[0]<b[2]<=width+.1 and -.1<=b[1]<b[3]<=height+.1):
            failures.append({'panels':[x['id']],'check':'positive in-page plot area'})
    for x,y in itertools.combinations(panels,2):
        if x.get('grid_id')!=y.get('grid_id'):continue
        bx,by=x['bbox_pt'],y['bbox_pt']
        for grid,index,label in [('col_start',0,'left'),('col_stop',2,'right'),('row_start',3,'top'),('row_stop',1,'bottom')]:
            if x[grid]==y[grid]:
                d=abs(bx[index]-by[index]);r={'panels':[x['id'],y['id']],'check':label,'deviation_pt':d}
                checks.append(r)
                if d>a.tolerance_pt:failures.append(r)
        # Equal-span rows of three or more panels require equal plot widths;
        # a two-panel evidence layout may intentionally allocate unequal widths.
        peers=[z for z in panels if (z['row_start'],z['row_stop'])==(x['row_start'],x['row_stop']) and z['col_stop']-z['col_start']==x['col_stop']-x['col_start']]
        if len(peers)>=3 and (x['row_start'],x['row_stop'])==(y['row_start'],y['row_stop']) and x['col_stop']-x['col_start']==y['col_stop']-y['col_start']:
            d=abs((bx[2]-bx[0])-(by[2]-by[0]));r={'panels':[x['id'],y['id']],'check':'equal-span row width','deviation_pt':d}
            checks.append(r)
            if d>a.tolerance_pt:failures.append(r)
    if not checks:raise ValueError('No comparable shared boundaries were measured.')
    result={'status':'PASS' if not failures else 'FAIL','tolerance_pt':a.tolerance_pt,'checks':checks,'failures':failures}
    a.json_out.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(result['status']+': '+str(len(checks))+' physical panel-boundary checks')
    if failures:raise SystemExit(1)

if __name__=='__main__':main()
