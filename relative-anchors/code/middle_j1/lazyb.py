"""A lazy copy of ylib.Builder: D^{-1} is computed per column sector, on demand (same entries as ylib.Builder)."""
from ylib import Builder, Mel, Nv, halfex, inv, wcol


class LazyBuilder(Builder):
    """D^{-1} computed per column sector, on demand."""

    def __init__(self, ring):
        super().__init__(ring)
        self._Dsec = {}
        self._bycol = {}
        for y in self.cols:
            self._bycol.setdefault(wcol(y), []).append(y)

    def Dinv_sec(self, key):
        if key not in self._Dsec:
            R_, k_, e = self.R, self.k, self.e
            ys = self._bycol[key]
            Bm = []
            for y in ys:
                row = []
                for yp in ys:
                    v = R_.zero
                    for tt in range(k_):
                        for tp in range(k_):
                            mv = Mel(e[tt], y, e[tp], yp)
                            if mv:
                                v += R_.one * mv * R_.mono(halfex(Nv(y), Nv(yp), e[tt], e[tp]))
                    row.append(v)
                Bm.append(row)
            Bi = inv(Bm, R_.zero, R_.one)
            self._Dsec[key] = {y: {yp: Bi[i][j] for j, yp in enumerate(ys) if Bi[i][j] != 0} for i, y in enumerate(ys)}
        return self._Dsec[key]

    def T(self, b, y):
        key = (b, y)
        if key not in self._T:
            acc = {}
            for z, v in self.MV(b, y).items():
                for zp, dv in self.Dinv_sec(wcol(z))[z].items():
                    acc[zp] = acc.get(zp, self.R.zero) + v * dv
            self._T[key] = acc
        return self._T[key]
