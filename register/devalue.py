import json,sys
def decode(arr):
    memo={}
    def h(i):
        if isinstance(i,int) and i<0: return None
        if i in memo: return memo[i]
        v=arr[i]
        if isinstance(v,dict):
            o={}; memo[i]=o
            for k,x in v.items(): o[k]=h(x)
            return o
        if isinstance(v,list):
            if v and isinstance(v[0],str) and v[0] in('Date','Set','Map','BigInt','RegExp'): return v[1] if len(v)>1 else None
            o=[]; memo[i]=o
            o.extend(h(x) for x in v); return o
        return v
    return h(0)
def load(p):
    d=json.load(open(p))
    return [decode(n['data']) for n in d['nodes'] if n and n.get('type')=='data']
