# Dictionary indexing (DI) of EBSD patterns from an Al-steel joint.
# Script to run on Idun computer cluster at NTNU.
#
# Håkon Wiik Ånes (hakon.w.anes@ntnu.no)
# 2022-05-28

import os

import dask
import matplotlib.pyplot as plt
import numpy as np
from orix import io, plot, sampling
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

# Sample orientation space
rot_dict = sampling.get_sample_fundamental(
    resolution=1.4,
    point_group=mp.phase.point_group,
    method="cubochoric",
)
print(rot_dict)

# Generate dictionary
sim_dict = mp.get_patterns(
    rotations=rot_dict,
    detector=detector,
    energy=energy,
    compute=False
)
print(sim_dict)

# Create a signal mask and check it
signal_mask = ~kp.filters.Window("circular", sig_shape).astype(bool)
p = s.inav[0, 0].data.compute()
fig, (ax0, ax1) = plt.subplots(ncols=2)
ax0.imshow(p * signal_mask, cmap="gray")
ax1.imshow(p * ~signal_mask, cmap="gray")
ax0.set_title("Not used")
ax1.set_title("Used")
fig.tight_layout()
fig.savefig(os.path.join(dir_data, "signal_mask.png"), **savefig_kwds)

# Index patterns
keep_n = 20
with dask.config.set(**{"array.slicing.split_large_chunks": False}):
    xmap = s.dictionary_indexing(
        dictionary=sim_dict,
        metric="ncc",
        keep_n=keep_n,
        signal_mask=signal_mask,
        rechunk=True,
    )
print(xmap)

# Save results to file
xmap.scan_unit = "um"
io.save(os.path.join(dir_data, f"xmap_{phase_name}.h5"), xmap, overwrite=True)

# Plot scores map
scores = xmap.scores[:, 0].reshape(xmap.shape)
plt.imsave(
    os.path.join(dir_data, f"maps_{phase_name}_ncc1.png"), scores, cmap="gray"
)

# Plot orientation similarity map
osm = kp.indexing.orientation_similarity_map(xmap, n_best=keep_n)
plt.imsave(
    os.path.join(dir_data, f"maps_{phase_name}_osm{keep_n}.png"),
    osm,
    cmap="gray"
)

# Both maps with colorbars
fig, (ax0, ax1) = plt.subplots(figsize=(12, 5), ncols=2)
im0 = ax0.imshow(scores, cmap="gray")
im1 = ax1.imshow(osm, cmap="gray")
fig.colorbar(im0, ax=ax0, label="Normalized cross correlation score")
fig.colorbar(im1, ax=ax1, label="Orientation similarity")
ax0.axis("off")
ax1.axis("off")
fig.tight_layout(w_pad=-1)
fig.savefig(
    os.path.join(dir_data, f"maps_{phase_name}_ncc1_osm{keep_n}_colorbar.png"),
    **savefig_kwds
)

# Histograms of both maps
bins_ncc = np.linspace(scores.min(), scores.max(), 100)
bins_osm = np.linspace(osm.min(), osm.max(), keep_n)
fig, (ax0, ax1) = plt.subplots(figsize=(9, 5), ncols=2)
ax0.hist(scores.ravel(), bins_ncc, color="C0", label="NCC from DI");
ax1.hist(osm.ravel(), bins_osm, color="C1", label="OSM");
ax0.set_xlim((0, 1))
ax1.set_xlim((0, keep_n))
ax0.set_xlabel("Normalized cross correlation score")
ax1.set_xlabel("Orientation similarity")
ax0.set_ylabel("Frequency")
ax1.set_ylabel("Frequency")
fig.tight_layout()
fig.savefig(
    os.path.join(dir_data, f"hist_{phase_name}_ncc1_osm{keep_n}.png"),
    **savefig_kwds
)

# Orientation maps
ckey = plot.IPFColorKeyTSL(xmap.phases[0].point_group)
ori = xmap.orientations
directions = Vector3d(((1, 0, 0), (0, 1, 0), (0, 0, 1)))

fig, axes = plt.subplots(figsize=(15, 5), ncols=3)
for ax, v, title in zip(axes, directions, ("x", "y", "z")):
    ckey.direction = v
    rgb = ckey.orientation2color(ori).reshape(xmap.shape + (3,))
    plt.imsave(
        os.path.join(dir_data, f"maps_{phase_name}_ipf{title}.png"), rgb
    )
    ax.imshow(rgb)
    ax.axis("off")
    ax.set_title(f"IPF {title}")
fig.tight_layout()
fig.savefig(
    os.path.join(dir_data, f"maps_{phase_name}_ipf.png"), **savefig_kwds
)
