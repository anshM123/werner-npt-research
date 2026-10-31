"""For NPT 3x3 filterings X of rho^2 decomposed as X = (a(x)b) rho (a(x)b)^* + S (S fitted by products): the spectrum and
   partial-transpose spectrum of the remainder R = X - (a(x)b) rho (a(x)b)^* (rank? PPT?), and the weight of each part."""
import numpy as np, torch
exec(open("check_kl1_structure.py", encoding="utf-8").read().split("def contraction_gap")[0])
torch.manual_seed(7); np.random.seed(7)
done = 0; seed = 0
while done < 3:
    seed += 1
    A = torch.randn(3, 9, dtype=torch.complex128, requires_grad=True); B = torch.randn(3, 9, dtype=torch.complex128, requires_grad=True)
    opt = torch.optim.Adam([A, B], lr=0.01)
    for it in range(20000):
        K = torch.kron(A, B); Y = K @ R2 @ K.conj().T; Y = (Y + Y.conj().T) / 2; Y = Y / torch.trace(Y).real
        f = torch.linalg.eigvalsh(pt(Y))[0]
        if f.item() < -0.01: break
        opt.zero_grad(); f.backward(); opt.step()
    if f.item() > -0.005: continue
    done += 1
    Yn = Y.detach().numpy()
    r, (a, b) = fit1(Yn, P=12, RS=10)
    W = np.kron(a, b) @ rho @ np.kron(a, b).conj().T
    Rm = Yn - W; Rm = (Rm + Rm.conj().T) / 2
    ev = np.linalg.eigvalsh(Rm); evp = np.linalg.eigvalsh(pt(torch.tensor(Rm)).numpy())
    print(f"target (PT-min {f.item():+.4f}): fit {r:.1e}; Tr parts: filtering {np.trace(W).real:.3f}, remainder {np.trace(Rm).real:.3f}")
    print(f"   remainder spectrum {np.round(ev / np.trace(Rm).real, 4)}")
    print(f"   remainder PT spectrum {np.round(evp / np.trace(Rm).real, 4)}", flush=True)
