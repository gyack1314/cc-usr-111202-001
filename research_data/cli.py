"""查看原始研究资料的已有摘要，并按截点推导受试者结局。"""
import argparse
import json
from .io import load
from .domain import summarize
from .endpoint import derive

def main():
    p=argparse.ArgumentParser();p.add_argument('command',choices=['inspect','derive']);p.add_argument('file');a=p.parse_args()
    value,digest=load(a.file)
    out={'source_sha256':digest}
    if a.command=='inspect':out['summary']=summarize(value)
    else:out['derivation']=derive(value)
    print(json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False))
if __name__=='__main__':main()
