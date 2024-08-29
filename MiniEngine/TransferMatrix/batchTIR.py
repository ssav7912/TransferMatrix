import numpy as np
from util import load_lut_from_file
from matplotlib import pyplot as plt

TIR_approx = np.fromfile("TIRapprox.bin").reshape((64,64,64))
TIR_naive = np.fromfile("naivereference.bin").reshape((64,64,64))
reference = load_lut_from_file("tm_TIR.bin", 3)

px = 1/plt.rcParams['figure.dpi']
text_kwargs = dict(ha='center', va='center', fontsize=24)
subplot_text = dict(fontsize=20)

for i in range(64):
    w, h = 64, 64
    print("Generating frame", i)
    naive_slice = plt.imsave(f"renders/naive/{i:2d}.png", TIR_naive[i], vmin=0.0, vmax=1.0, cmap="Reds_r")
    approx_slice = plt.imsave(f"renders/approx/{i:2d}.png", TIR_approx[i], vmin=0.0, vmax=1.0, cmap="Reds_r")
    ref_slice = plt.imsave(f"renders/reference/{i:2d}.png", reference[i], vmin=0.0, vmax=1.0, cmap="Reds_r")

    ext=[0,1,0,1]
    f, axarr = plt.subplots(3, 1, sharex=True, sharey=True, figsize=(1080*px, 1920*px))
    f.suptitle("Approximate Analytic Form For Total Internal Reflection", **text_kwargs, y=0.95)
    nplot = axarr[1].imshow(TIR_naive[i], vmin=0.0, vmax=1.0, extent=ext, cmap="Reds_r", )
    approxplot = axarr[2].imshow(TIR_approx[i], vmin=0.0, vmax=1.0, extent=ext, cmap="Reds_r")
    refplot = axarr[0].imshow(reference[i], vmin=0.0, vmax=1.0, extent=ext, cmap="Reds_r")

    axarr[1].set_title("Naive TIR Estimator (Belcour, 2018)", **subplot_text)
    axarr[0].set_title("Reference TIR Estimator (Belcour, 2018)", **subplot_text)
    axarr[2].set_title("Approximate Analytical TIR", **subplot_text)
    axarr[2].set_xlabel("Interface Alpha (Roughness)", **subplot_text)

    for subplot in axarr:
        #subplot.set_xlabel("Interface Alpha (Roughness)", **subplot_text)
        subplot.set_ylabel("Ratio between layer IORs", **subplot_text)
        subplot.xaxis.set_tick_params(labelsize='x-large')
        subplot.yaxis.set_tick_params(labelsize='x-large')
    
    f.text(x=0.1, y=0.92, s=f"Incident Angle {i/64:.2f}", ha='left', va='center', fontsize=24)

    col = f.colorbar(nplot, ax=axarr, orientation='vertical', fraction=0.1)
    col.set_label(label="Proportion of Reflected Energy", **subplot_text)
    col.ax.tick_params(labelsize='x-large')
    

    f.savefig(f"renders/graph{i:002d}.png")
    plt.close(f)