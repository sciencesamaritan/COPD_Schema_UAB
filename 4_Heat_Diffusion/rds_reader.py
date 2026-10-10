"""Minimal reader for R serialization (RDS, XDR v2/v3): lists, atomic vectors, attributes, ALTREP compact seqs."""
import gzip, struct, numpy as np

class RObj:
    def __init__(self, value, attrs=None, rtype=None): self.value, self.attrs, self.rtype = value, attrs or {}, rtype
    def __repr__(self): return f"RObj({self.rtype}, attrs={list(self.attrs)})"

class Reader:
    def __init__(self, b): self.b, self.p, self.refs = b, 0, []
    def i(self):
        v = struct.unpack(">i", self.b[self.p:self.p + 4])[0]; self.p += 4; return v
    def n(self):
        l = self.i()
        if l == -1: hi, lo = self.i(), self.i(); l = (hi << 32) + lo
        return l
    def charsxp(self, flags=None):
        if flags is None: flags = self.i()
        l = self.i()
        if l == -1: return None
        s = self.b[self.p:self.p + l].decode("utf-8", "replace"); self.p += l; return s
    def attrs_pairlist(self):
        obj = self.item(); return obj if isinstance(obj, dict) else {}
    def item(self):
        flags = self.i(); t = flags & 0xFF; has_attr = flags & (1 << 9); has_tag = flags & (1 << 10)
        if t == 254: return None                                  # NILVALUE
        if t in (253, 242, 251, 250, 249, 248, 247, 246): return None  # env/missing/unbound markers
        if t == 255:                                               # REFSXP
            idx = flags >> 8; idx = idx if idx else self.i(); return self.refs[idx - 1]
        if t == 1:                                                 # SYMSXP
            s = self.item(); self.refs.append(s); return s
        if t == 9: return self.charsxp(flags)
        if t in (2, 6, 17):                                        # LISTSXP (pairlist) -> dict
            out = {}
            while True:
                if has_attr: self.item()
                tag = self.item() if has_tag else None
                car = self.item(); out[tag if tag is not None else len(out)] = car
                flags = self.i(); t = flags & 0xFF; has_attr = flags & (1 << 9); has_tag = flags & (1 << 10)
                if t == 254: break
                if t not in (2, 6, 17): raise ValueError(f"bad pairlist cdr type {t}")
            return out
        if t == 4:                                                 # ENVSXP
            self.i(); env = {}; self.refs.append(env)
            self.item(); self.item(); self.item(); self.item(); return env
        if t == 238:                                               # ALTREP
            info = self.item(); state = self.item(); attr = self.item()
            cls = info.get(0) if isinstance(info, dict) else None
            if cls == "compact_intseq":
                v = state.value if isinstance(state, RObj) else state; n, a, s = int(v[0]), int(v[1]), int(v[2])
                return RObj(np.arange(a, a + s * n, s), attr if isinstance(attr, dict) else {}, "int")
            if cls == "compact_realseq":
                v = state.value; n, a, s = int(v[0]), v[1], v[2]; return RObj(a + s * np.arange(n), {}, "real")
            if isinstance(state, RObj) and state.rtype == "list": return state.value[0]
            return state
        if t in (10, 13):
            l = self.n(); v = np.frombuffer(self.b, ">i4", l, self.p).astype(np.int64); self.p += 4 * l
            return RObj(v, self.attrs_pairlist() if has_attr else {}, "int" if t == 13 else "lgl")
        if t == 14:
            l = self.n(); v = np.frombuffer(self.b, ">f8", l, self.p).astype(float); self.p += 8 * l
            return RObj(v, self.attrs_pairlist() if has_attr else {}, "real")
        if t == 16:
            l = self.n(); v = [self.charsxp() for _ in range(l)]
            return RObj(v, self.attrs_pairlist() if has_attr else {}, "str")
        if t in (19, 20):
            l = self.n(); v = [self.item() for _ in range(l)]
            return RObj(v, self.attrs_pairlist() if has_attr else {}, "list")
        raise ValueError(f"unsupported SEXP type {t} at {self.p}")

def read_rds(path):
    b = gzip.open(path).read() if open(path, "rb").read(2) == b"\x1f\x8b" else open(path, "rb").read()
    r = Reader(b); fmt = b[:2]; assert fmt == b"X\n", fmt; r.p = 2
    ver = r.i(); r.i(); r.i()
    if ver == 3: nl = r.i(); r.p += nl
    return r.item()

def names(o): return o.attrs.get("names").value if isinstance(o, RObj) and "names" in o.attrs else None
