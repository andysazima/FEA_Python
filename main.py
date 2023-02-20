# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 19:07:40 2022
@author: andys

Main run file
"""
from fem_solver import FEMSolver
from sphere import Sphere
import time
import matplotlib.pyplot as plt
plt.rcParams["figure.dpi"] = 500

start = time.time()


sphere = Sphere()
solver = FEMSolver(sphere)
solver.spherical_pressure()

end = time.time()

print(f"Time Elapsed = {end - start:.4g} sec")

prob_text = rf"E = {sphere.E}" + "\n" + \
            rf"$\nu$ = {sphere.nu}" + "\n" + \
            rf"$\rho$ = {sphere.rho}" + "\n" + \
            rf"$R_i$ = {sphere.R_i}" + "\n" + \
            rf"$R_o$ = {sphere.R_o}" + "\n" + \
            rf"# Ele = {sphere.n_ele}" + "\n" + \
            rf"$\alpha_0$ = {sphere.alpha_0}" + "\n" + \
            rf"$\alpha_1$ = {sphere.alpha_1}" + "\n" + \
            rf"$P_0$ = {sphere.P_0}" + "\n\n" + \
            rf"$\beta$ = {solver.beta}" + "\n" + \
            rf"$\gamma$ = {solver.gamma}" + "\n" + \
            r"$t_{tot}$ = " + f"{solver.t_tot}" + "\n" + \
            r"$\Delta t$ = " + f"{solver.dt:.4g}"

fig, ax = plt.subplots(3, 1)
fig.text(0.92, 0.35, prob_text, fontsize=11)
ax[0].set_title("Time Histories")
fig.set(figwidth=13,
        figheight=8)
ax[0].axhline(y=0, color='k')
ax[0].plot(solver.t_hist, sphere.d[0, :],  label=r"$r = R_i$")
ax[0].plot(solver.t_hist, sphere.d[int(
    sphere.n_ele/2), :], label=r"$r = R_o/2$")
ax[0].plot(solver.t_hist, sphere.d[-1, :], label=r"$r = R_o$")
ax[0].set(ylabel="Displacement [L]")
ax[0].margins(x=0)
ax[0].grid(True)

ax[1].axhline(y=0, color='k')
ax[1].plot(solver.t_hist, sphere.sig[0, :, 0],  label=r"$r = R_i$")
ax[1].plot(solver.t_hist, sphere.sig[0, :, int(
    sphere.n_ele/2)], label=r"$r = R_o/2$")
ax[1].plot(solver.t_hist, sphere.sig[0, :, -1], label=r"$r = R_o$")
ax[1].set(ylabel=r"Radial Stress [F/L$^2$]")
ax[1].margins(x=0)
ax[1].grid(True)

ax[2].axhline(y=0, color='k')
ax[2].plot(solver.t_hist, sphere.sig[1, :, 0],  label=r"$r = R_i$")
ax[2].plot(solver.t_hist, sphere.sig[1, :, int(
    sphere.n_ele/2)], label=r"$r = R_o/2$")
ax[2].plot(solver.t_hist, sphere.sig[1, :, -1], label=r"$r = R_o$")
ax[2].set(xlabel="Time [T]",
          ylabel="Circumferential Stress [F/L$^2$]")
ax[2].margins(x=0)
ax[2].grid(True)
ax[2].legend(bbox_to_anchor=(0.5, -0.25), loc='upper center', ncol=3)

filename = "Time Histories t_" + str(int(sphere.R_o - sphere.R_i))
fig.savefig(filename, bbox_inches="tight")


# # fig, ax = plt.subplots(3,1)
# fig, ax = plt.subplots(1,1)
# fig.text(0.92, 0., prob_text, fontsize=11)
# # ax.set_title("Time Histories")
# fig.set(figwidth  = 5,
#         figheight = 2.5)
# ax.axhline(y = 0, color = 'k')
# ax.plot(solver.t_hist, sphere.sig[1,:,0,0],  label = r"$r = R_i$")
# ax.plot(solver.t_hist, sphere.sig[1,:,int(sphere.n_ele/2),0],
# label = r"$r = R_o/2$")
# ax.plot(solver.t_hist, sphere.sig[1,:,-1,-1], label = r"$r = R_o$")
# ax.set(xlabel = "Time [T]",
#         ylabel = "Circumferential Stress [F/L$^2$]")
# ax.margins(x = 0)
# ax.grid(True)
# ax.legend(bbox_to_anchor = (0.5,-0.2), loc='upper center', ncol = 3)

# filename = "Thickness Study t_" + str(int(sphere.R_o - sphere.R_i))
# fig.savefig(filename, bbox_inches="tight")


# # fig, ax = plt.subplots(3,1)
# r_plot = np.reshape(sphere.r_xi[::2], (sphere.n_ele, 1))
# fig, ax = plt.subplots(1,1)
# # fig.text(0.92, 0., prob_text, fontsize=11)
# # ax.set_title("Time Histories")
# fig.set(figwidth  = 5,
#         figheight = 2.5)
# ax.axhline(y = 0, color = 'k')
# ax.plot(r_plot, sphere.sig[1,-1,:],  label = r"$r = R_i$")
# ax.set(xlabel = "Radial Position Along Cross Section [L]",
#         ylabel = "Circumferential Stress [F/L$^2$]")
# ax.margins(x = 0.005)
# ax.grid(True)

# filename = "Thickness Study Stress Dist t_" + \
# str(int(sphere.R_o - sphere.R_i))
# fig.savefig(filename, bbox_inches="tight")


# print(str(sphere.n_ele) + ",")
# print(str(np.max(sphere.d)) + ",")

# Convergence = [1,
# 0.11379566955374931,
# 2,
# 0.14211780535712026,
# 3,
# 0.1463442679195739,
# 4,
# 0.1493950258378489,
# 6,
# 0.15183169759723925,
# 10,
# 0.1490328530820643,
# 20,
# 0.15120651696921517,
# 30,
# 0.1538437103241082,
# 40,
# 0.15531531011968813,
# 60,
# 0.1541041890291717,
# 100,
# 0.15308687092936957,
# 200,
# 0.15345810076936411,
# 300,
# 0.15332957005082593,
# 400,
# 0.15332978732173436,
# 500,
# 0.15332373237992447,
# 1000,
# 0.15334622718679608,
# ]
# conv_plot = np.reshape(Convergence, (2, int(len(Convergence)/2)),
# order = "F")

# fig, ax = plt.subplots(1,1)
# fig.text(0.92, 0., prob_text, fontsize=11)
# # ax.set_title("Time Histories")
# fig.set(figwidth  = 5,
#         figheight = 2.5)
# ax.grid(True)
# ax.axhline(y = 0, color = 'k')
# ax.scatter(conv_plot[0,:], conv_plot[1,:])
# ax.set(xlabel = "Number of Elements [-]",
#        ylabel = "Maximum Displacement [L]",
#        xscale = "log")


# filename = "Num Ele" + str(int(sphere.R_o - sphere.R_i))
# fig.savefig(filename, bbox_inches="tight")
