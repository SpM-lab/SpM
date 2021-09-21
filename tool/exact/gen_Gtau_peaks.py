import math
import numpy as np
import scipy.integrate


def gaussian(x, x0, sigma):
    return 1.0 / (math.sqrt(math.pi) * sigma) * math.exp(-(((x - x0) / sigma) ** 2))


def lorentzian(x, x0, gamma):
    return (1.0 / math.pi) * gamma / ((x - x0) ** 2 + gamma ** 2)


def semicircle(x, x0, width):
    if abs(x - x0) >= width:
        return 0
    return (2.0 / math.pi / width ** 2) * (width ** 2 - (x - x0) ** 2) ** 0.5


class Spectrum:
    func = {"G": gaussian, "L": lorentzian, "S": semicircle}

    def __init__(self, stat="fermion"):
        self.stat = stat
        self.peaks = []

    def add_peak(self, peak_type, position, width, intensity):
        self.peaks.append([peak_type, position, width, intensity])
        if self.stat == "boson":
            self.peaks.append([peak_type, -position, width, -intensity])

    def dos(self, w):
        s = 0
        for which, pos, width, A in self.peaks:
            s += A * Spectrum.func[which](w, pos, width)
        return s

    def save(self, filename="Gtau.in.dos", omega=[0.0]):
        with open(filename, "w") as f:
            for which, pos, width, A in self.peaks:
                print(
                    f"# (type, pos, width, intensity) = ({which}, {pos}, {width}, {A})",
                    file=f,
                )
            print("# Nomega = ", len(omega), file=f)
            print("# omega_max = ", omega[-1], file=f)
            for w in omega:
                print(w, self.dos(w), file=f)


class GreenFunction:
    def __init__(self, beta, Ntau):
        self.set_param(beta, Ntau)

    def set_param(self, beta, Ntau):
        self.beta = beta
        self.Ntau = Ntau
        self.tau = np.linspace(0, beta, Ntau)
        self.Gtau = np.zeros(self.Ntau)
        self.calculated = False

    def calc_Gtau(self, spectrum, omega_min, omega_max, stat="fermion"):
        if stat == "fermion":
            self.kernel = self.kernel_f
        elif stat == "boson":
            self.kernel = self.kernel_b
        else:
            raise RuntimeError(
                f"Unknown statistics type: {stat} (fermion or boson is available)"
            )

        def integrand(w, t):
            r = spectrum.dos(w)
            if r == 0.0:
                return 0.0
            else:
                return r * self.kernel(t, w)

        for i, t in enumerate(self.tau):
            self.Gtau[i] = scipy.integrate.quad(
                integrand, omega_min, omega_max, args=t
            )[0]
        if stat == "fermion":
            self.Gtau *= -1.0
        self.calculated = True
        self.omega_max = omega_max
        self.omega_min = omega_min

    def kernel_f(self, tau, omega):
        if omega >= 0:
            return math.exp(-tau * omega) / (1.0 + math.exp(-self.beta * omega))
        else:
            return math.exp((self.beta - tau) * omega) / (
                1.0 + math.exp(self.beta * omega)
            )

    def kernel_b(self, tau, omega):
        bw = self.beta * omega
        if omega > 0:
            return -np.exp(-tau * omega) / np.expm1(-bw)
        else:
            return np.exp((self.beta - tau) * omega) / np.expm1(bw)

    def save(self, filename="Gtau.in", noise=0.0):
        if not self.calculated:
            print("Gtau is not calculated")
            return

        noised = np.random.normal(loc=self.Gtau, scale=noise)

        with open(filename, "w") as f:
            print(f"#beta={self.beta}", file=f)
            print(f"#omega_max={self.omega_max}", file=f)
            print(f"#omega_min={self.omega_min}", file=f)
            print(f"#noise_sigma={noise}", file=f)
            print(f"#tau/beta, w/ noise, w/o noise", file=f)
            for t, Gn, G in zip(self.tau, noised, self.Gtau):
                print(t / self.beta, Gn, G, file=f)


if __name__ == "__main__":
    import os
    import argparse
    import toml

    parser = argparse.ArgumentParser(
        description="To generate Gtau with white noise from a spectrum", add_help=True
    )

    parser.add_argument("input", type=argparse.FileType("r"), help="Input TOML file")

    parser.add_argument(
        "-n",
        dest="nsamples",
        default=1,
        type=int,
        help="the number of samples to be generated (default=%(default)s)",
    )

    parser.add_argument("-q", action="store_true", help="do quietly")

    args = parser.parse_args()

    VERBOSE = not args.q

    param = toml.load(args.input)

    if VERBOSE:
        print("parameter file is successfully loaded")

    outdir = param["parameter"].get("outputdir", "work")
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        if VERBOSE:
            print(f"make output directory ({outdir})")

    if not os.path.isdir(outdir):
        raise RuntimeError(f"{outdir} is not a directory.")

    nsamples = args.nsamples

    beta = param["parameter"]["beta"]
    Ntau = param["parameter"]["Ntau"]
    omega_min = param["parameter"]["omega_min"]
    omega_max = param["parameter"]["omega_max"]
    Nomega = param["parameter"]["Nomega"]
    noise = param["parameter"]["noise"]
    stat = param["parameter"].get("stat", "fermion")

    if VERBOSE:
        print(f"stat = {stat}")
        print(f"beta = {beta}")
        print(f"Ntau = {Ntau}")
        print(f"Nomega = {Nomega}")
        print(f"omega_min = {omega_min}")
        print(f"omega_max = {omega_max}")
        print(f"noise = {noise}")

    seed_base = param["parameter"].get("seed_base", 314159)
    seed_step = param["parameter"].get("seed_step", 11)

    spectrum = Spectrum(stat)
    for peak in param["peaks"]:
        spectrum.add_peak(**peak)

    omega = np.linspace(omega_min, omega_max, num=Nomega)
    dosfilename = os.path.join(outdir, "Gtau.in.dos")
    spectrum.save(dosfilename, omega)
    if VERBOSE:
        print(f"Spectrum is saved ({dosfilename})")

    gf = GreenFunction(beta, Ntau)
    gf.calc_Gtau(spectrum, omega_min, omega_max, stat=stat)
    if VERBOSE:
        print("Exact Green's function is calculated")

    if nsamples == 1:
        gtaufilename = os.path.join(outdir, "Gtau.in")
        np.random.seed(seed_base)
        gf.save(gtaufilename, noise)
        if VERBOSE:
            print(f"Green's function is saved ({gtaufilename})")
    else:
        for i in range(nsamples):
            gtaufilename = os.path.join(outdir, f"Gtau_{i}.in")
            np.random.seed(seed_base + i * seed_step)
            gf.save(gtaufilename, noise)
            if VERBOSE:
                print(f"{i+1}/{nsamples} Green's function is saved ({gtaufilename})")
