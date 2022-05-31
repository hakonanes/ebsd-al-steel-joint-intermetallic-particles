# Orientation refinement of EBSD patterns from an Al-steel joint.
# Script to run on Idun computer cluster at NTNU.
#
# Håkon Wiik Ånes (hakon.w.anes@ntnu.no)
# 2022-05-28

import os

import matplotlib.pyplot as plt
import numpy as np
from orix import io, plot
from orix.vector import Vector3d

import kikuchipy as kp


# Directories
# Dataset naming (a-c) = (I-III)
dset = "a"
phase_name = "ferrite"
dir_mp = "/cluster/work/hakonwii/emdata/crystal_data"
dir_data = os.path.join("/cluster/home/hakonwii/data/tina", dset, phase_name)

# Projection center (PC)
pcs = {
    "a": (0.42586528, 0.19035024, 0.51163327),
    "b": (0.42865885, 0.18983122, 0.51430069),
    "c": (0.42287091, 0.18307211, 0.50742970),
    "d": (0.42192298, 0.18911612, 0.51270242),
}

# Matplotlib parameters
savefig_kwds = dict(dpi=200, bbox_inches="tight", pad_inches=0)

# Load data lazily
s = kp.load(os.path.join(dir_data, f"pattern_sda_{phase_name}.h5"), lazy=True)
sig_shape = s.axes_manager.signal_shape[::-1]
print(s)

# Detector's view of the sample
detector = kp.detectors.EBSDDetector(
    shape=sig_shape,
    pc=pcs[dset],
    sample_tilt=70,
    convention="bruker",
)
print(detector)

# Check projection center
fig, _ = detector.plot(pattern=s.inav[0, 0].data.compute(), return_fig_ax=True)
fig.savefig(os.path.join(dir_data, "detector.png"), **savefig_kwds)

# Load master pattern
energy = 20
mp = kp.load(
    os.path.join(dir_mp, phase_name, f"{phase_name}_mc_mp_20kv.h5"),
    projection="lambert",
    energy=energy,
    hemisphere="both",
)
mp.phase.name = phase_name
print(mp)

# Load initial orientation results
xmap = io.load(os.path.join(dir_data, f"xmap_{phase_name}.h5"))
print(xmap)

# Refinement
xmap_ref = s.refine_orientation(
    xmap=xmap,
    detector=detector,
    master_pattern=mp,
    energy=energy,
    compute=True,
)
print(xmap_ref)
io.save(os.path.join(dir_data, f"xmap_ref_{phase_name}.h5"), xmap_ref)

# Plot normalized cross-correlation scores
scores1d_di = xmap.scores[:, 0]
scores2d_di = scores1d_di.reshape(xmap.shape)

scores1d = xmap_ref.scores
scores2d = scores1d.reshape(xmap_ref.shape)
plt.imsave(
    os.path.join(dir_data, f"maps_{phase_name}_ncc_ref.png"),
    arr=scores2d,
    cmap="gray"
)

ncc_di_min = np.min(scores1d)
ncc_di_max = np.max(scores1d)
ncc_ori_ref_min = np.min(scores1d_di)
ncc_ori_ref_max = np.max(scores1d_di)
vmin = min([ncc_di_min, ncc_ori_ref_min])
vmax = max([ncc_di_max, ncc_ori_ref_max])

ncc_after_label = f"NCC after ref. ori."

fig, ax = plt.subplots(ncols=2, figsize=(10, 3))
im0 = ax[0].imshow(scores2d_di, vmin=vmin, vmax=vmax)
im1 = ax[1].imshow(scores2d, vmin=vmin, vmax=vmax)
fig.colorbar(im0, ax=ax[0], label="NCC from DI")
fig.colorbar(im1, ax=ax[1], label=ncc_after_label)
for a in ax:
    a.axis("off")
fig.tight_layout(w_pad=-1)
fig.savefig(
    os.path.join(dir_data, f"maps_{phase_name}_ncc_comparison_ref.png"),
    **savefig_kwds
)

bins = np.linspace(vmin, vmax, 100)

fig, ax = plt.subplots(figsize=(9, 5))
ax.hist(scores1d_di, bins, alpha=0.5, label="NCC from DI");
ax.hist(scores1d, bins, alpha=0.5, label=ncc_after_label);
ax.set_xlabel("Normalized cross correlation (NCC) scores")
ax.set_ylabel("Frequency")
ax.legend()
fig.tight_layout()
fig.savefig(
    os.path.join(dir_data, f"hist_{phase_name}_ncc_comparison_ref.png"),
    **savefig_kwds
)

# Plot maps of IPFZ and disorientation angle
ori = xmap.orientations
ori_ref = xmap_ref.orientations

ckey = plot.IPFColorKeyTSL(xmap.phases[0].point_group)
ckey.direction = Vector3d.xvector()
rgb_z = ckey.orientation2color(ori)
rgb_z_ref = ckey.orientation2color(ori_ref)

mori_angle = np.rad2deg(ori.angle_with(ori_ref).data)
vmin, vmax = np.percentile(mori_angle, q=(1, 99))

fig, ax = plt.subplots(ncols=3, figsize=(14, 3))
ax[0].imshow(rgb_z.reshape(xmap.shape + (3,)))
ax[1].imshow(rgb_z_ref.reshape(xmap.shape + (3,)))
im2 = ax[2].imshow(mori_angle.reshape(xmap.shape), vmin=vmin, vmax=vmax)
fig.colorbar(im2, ax=ax[2], label=r"Disorientation angle $\omega$ [$^{\circ}$]")
for a in ax:
    a.axis("off")
fig.tight_layout(w_pad=-15)
fig.savefig(
    os.path.join(dir_data, f"maps_{phase_name}_ipfx_comparison_ref.png"),
    **savefig_kwds
)

# Plot histogram of misorientation angle between DI and refined orientations
fig, ax = plt.subplots(figsize=(9, 5))
ax.hist(mori_angle, bins=50, range=(vmin, vmax))
ax.set_xlim((vmin, vmax))
ax.set_xlabel("Misorientation angle DI -> ref. $\omega$ [$^{\circ}$]")
ax.set_ylabel("Frequency")
fig.tight_layout()
fig.savefig(
    os.path.join(dir_data, f"mori_angle_di_ref_{phase_name}.png"),
    **savefig_kwds
)

# Plot orientation maps
directions = Vector3d(((1, 0, 0), (0, 1, 0), (0, 0, 1)))
map_shape = xmap_ref.shape

fig, axes = plt.subplots(figsize=(15, 5), ncols=3)
for ax, v, title in zip(axes, directions, ("x", "y", "z")):
    ckey.direction = v
    rgb = ckey.orientation2color(ori_ref).reshape(map_shape + (3,))
    plt.imsave(
        os.path.join(dir_data, f"maps_{phase_name}_ipf{title}_ref.png"), rgb
    )
    ax.imshow(rgb)
    ax.axis("off")
    ax.set_title(f"IPF {title}")
fig.tight_layout()
fig.savefig(
    os.path.join(dir_data, f"maps_{phase_name}_ipf_ref.png"), **savefig_kwds
)
