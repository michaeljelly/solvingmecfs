import requests, io
class HTTPFile(io.RawIOBase):
    def __init__(self,url,block=1<<20):
        self.url=url; self.s=requests.Session(); self.pos=0; self.block=block; self.cache={}
        r=self.s.head(url,allow_redirects=True); self.size=int(r.headers['Content-Length'])
    def seekable(self): return True
    def readable(self): return True
    def tell(self): return self.pos
    def seek(self,off,whence=0):
        self.pos = off if whence==0 else (self.pos+off if whence==1 else self.size+off); return self.pos
    def _blk(self,i):
        if i not in self.cache:
            if len(self.cache)>64: self.cache.pop(next(iter(self.cache)))
            a=i*self.block; b=min(self.size,a+self.block)-1
            for t in range(5):
                try:
                    r=self.s.get(self.url,headers={'Range':f'bytes={a}-{b}'},timeout=120); r.raise_for_status(); self.cache[i]=r.content; break
                except Exception as e:
                    if t==4: raise
        return self.cache[i]
    def read(self,n=-1):
        if n is None or n<0: n=self.size-self.pos
        n=min(n,self.size-self.pos); out=bytearray()
        while n>0:
            i=self.pos//self.block; o=self.pos%self.block; d=self._blk(i)[o:o+n]
            out+=d; self.pos+=len(d); n-=len(d)
        return bytes(out)
    def readinto(self,b):
        d=self.read(len(b)); b[:len(d)]=d; return len(d)
