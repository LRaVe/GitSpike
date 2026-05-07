# -*- coding: utf-8 -*-
"""
Created on Wed May  6 20:53:14 2026

@author: maxim
"""


# -----------------------------
# SPIKE trains (same as MATLAB)
# -----------------------------
spikes = [
    [0, 1, 2, 4, 7, 10],
    [0, 3, 4, 6, 10],
    [0, 1.9, 3.9, 7, 10]
]

# -----------------------------
# Time window
# -----------------------------
t_min = 0
t_max = 10



    
import numpy as np


def spike_dist_2x2(spikes1, spikes2, t_min, t_max):

    spikes1 = np.array(spikes1, dtype=float)
    spikes2 = np.array(spikes2, dtype=float)

    profile = []
    profile.append([t_min, 0.0])

    t_all = np.sort(np.concatenate([spikes1, spikes2]))

    for k in range(2, len(t_all)-2):

        t = t_all[k]

        if np.sum(t_all == t) > 1:
            profile.append([t, 0.0])
            continue

        # -------------------------
        # CASE spikes1 event
        # -------------------------
        if np.any(spikes1 == t):

            idx1 = np.where(spikes1 == t)[0][0]
            idx2 = np.where(spikes2 <= t)[0][-1]

            ISI2 = spikes2[idx2+1] - spikes2[idx2]

            dtp2 = np.min(np.abs(spikes2[idx2] - spikes1))
            dtf2 = np.min(np.abs(spikes2[idx2+1] - spikes1))

            xp2 = t - spikes2[idx2]
            xf2 = spikes2[idx2+1] - t

            S2 = (dtp2*xf2 + dtf2*xp2) / ISI2

            # LEFT
            ISI1_L = spikes1[idx1] - spikes1[idx1-1]
            dtf1 = np.min(np.abs(spikes1[idx1] - spikes2))
            xp1 = t - spikes1[idx1-1]

            S1 = (dtf1*xp1) / ISI1_L

            left = ((S1*ISI2) + (S2*ISI1_L)) / (2*(np.mean([ISI1_L, ISI2])**2))

            # RIGHT
            ISI1_R = spikes1[idx1+1] - spikes1[idx1]
            dtp1 = np.min(np.abs(spikes1[idx1] - spikes2))
            xf1 = spikes1[idx1+1] - t

            S1 = (dtp1*xf1) / ISI1_R

            right = ((S1*ISI2) + (S2*ISI1_R)) / (2*(np.mean([ISI1_R, ISI2])**2))

            profile.append([t, left])
            profile.append([t, right])

        # -------------------------
        # CASE spikes2 event
        # -------------------------
        else:

            idx1 = np.where(spikes1 <= t)[0][-1]
            idx2 = np.where(spikes2 == t)[0][0]

            ISI1 = spikes1[idx1+1] - spikes1[idx1]

            dtp1 = np.min(np.abs(spikes1[idx1] - spikes2))
            dtf1 = np.min(np.abs(spikes1[idx1+1] - spikes2))

            xp1 = t - spikes1[idx1]
            xf1 = spikes1[idx1+1] - t

            S1 = (dtp1*xf1 + dtf1*xp1) / ISI1

            # LEFT
            ISI2_L = spikes2[idx2] - spikes2[idx2-1]
            dtf2 = np.min(np.abs(spikes2[idx2] - spikes1))
            xp2 = t - spikes2[idx2-1]

            S2 = (dtf2*xp2) / ISI2_L

            left = ((S1*ISI2_L) + (S2*ISI1)) / (2*(np.mean([ISI1, ISI2_L])**2))

            # RIGHT
            ISI2_R = spikes2[idx2+1] - spikes2[idx2]
            dtp2 = np.min(np.abs(spikes2[idx2] - spikes1))
            xf2 = spikes2[idx2+1] - t

            S2 = (dtp2*xf2) / ISI2_R

            right = ((S1*ISI2_R) + (S2*ISI1)) / (2*(np.mean([ISI1, ISI2_R])**2))

            profile.append([t, left])
            profile.append([t, right])

    profile.append([t_max, 0.0])

    profile = np.array(profile, dtype=float)

    t = profile[:, 0]
    S = profile[:, 1]

    D = np.trapz(S, t) / (t_max - t_min)

    return D, profile


def spike_dist_N(spikes, t_min, t_max):

    N = len(spikes)
    D_matrix = np.zeros((N, N))
    profiles = []

    # pairwise computation
    for i in range(N):
        for j in range(i+1, N):

            d, prof = spike_dist_2x2(spikes[i], spikes[j], t_min, t_max)

            D_matrix[i, j] = d
            D_matrix[j, i] = d

            profiles.append(prof)

    D_global = np.mean(D_matrix[np.triu_indices(N, 1)])

    # global time axis
    t_all = np.unique(np.concatenate([p[:, 0] for p in profiles]))
    t_all.sort()

    profile_global = []

    for t in t_all:

        left_vals = []
        right_vals = []
        single_vals = []
        has_disc = False

        for P in profiles:

            idx = np.where(P[:, 0] == t)[0]

            if len(idx) == 2:
                has_disc = True
                left_vals.append(P[idx[0], 1])
                right_vals.append(P[idx[1], 1])

            elif len(idx) == 1:
                single_vals.append(P[idx[0], 1])

            else:
                before = np.where(P[:, 0] < t)[0]
                after = np.where(P[:, 0] > t)[0]

                if len(before) > 0 and len(after) > 0:
                    i1 = before[-1]
                    i2 = after[0]

                    t1, t2 = P[i1, 0], P[i2, 0]
                    s1, s2 = P[i1, 1], P[i2, 1]

                    interp = s1 + (s2 - s1) * (t - t1) / (t2 - t1)
                    single_vals.append(interp)

        if has_disc:

            if single_vals:
                left_vals += single_vals
                right_vals += single_vals

            if left_vals:
                profile_global.append([t, np.mean(left_vals)])
            if right_vals:
                profile_global.append([t, np.mean(right_vals)])

        else:
            if single_vals:
                profile_global.append([t, np.mean(single_vals)])

    profile_global = np.array(profile_global, dtype=float)

    # plot
    import matplotlib.pyplot as plt

    plt.figure()
    plt.plot(profile_global[:, 0], profile_global[:, 1])
    plt.xlim(t_min, t_max)
    plt.ylim(0, 1)
    plt.title(f"SPIKE-distance = {D_global:.4f}")
    plt.xlabel("Time")
    plt.ylabel("SPIKE distance")
    plt.show()

    # matrix
    plt.figure()
    plt.imshow(D_matrix, cmap="jet", interpolation="nearest")
    plt.colorbar()
    plt.xticks(np.arange(N), np.arange(1, N+1))
    plt.yticks(np.arange(N), np.arange(1, N+1))
    plt.title(f"SPIKE-distance = {D_global:.4f}")
    plt.xlabel("SPIKE Trains")
    plt.ylabel("SPIKE Trains")
    plt.show()

    return D_global, profile_global, D_matrix





# -----------------------------
# Run SPIKE-distance
# -----------------------------
D, profile, M = spike_dist_N(spikes, t_min, t_max)

# -----------------------------
# Display results
# -----------------------------
print("\nGlobal SPIKE-distance:", D)
print("\nMatrix D:")
print(M)








