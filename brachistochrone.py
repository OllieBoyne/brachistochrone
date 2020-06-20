from matplotlib import pyplot as plt
from viseng.figure import FuncFig
import numpy as np
from time import perf_counter
from tqdm import tqdm
import sys
sys.path.append("../")

from brachistochrone_analytical import *

## To produce a brachistochrone, the time taken for a ball to get from one end to the other must be minimised.
# curve is y(x), where
# y(0) = H
# y(L) = 0
# t = integral between 0 and L of sqrt( (y'' + 1) / 2gy) dx


font = "Arial"

hyperparams = {
	"newton": dict(method = "newton", nit = 200, delta = 0.0001),
	"grad_descent": dict(method="grad_descent", nit = 1000, delta = 0.001, base_lr=0.05
						 )
}

# brachisochrone params
H = 1
L = 2

base = 5
log_base = lambda n, base: np.log(n) / np.log(base)
n_samples = 25
X = np.logspace(-2, log_base(L, base), n_samples-1, base = base) # sample x in a log scale. This manages the fact that more time is spent earlier on.

# sample geometrically instead
# X = np.cumsum(np.arange(1, n_samples+1))
# X = (X * L / X.max())[:-1] # normalise so 0 < X < L

X = np.concatenate([ [0], X])

def compute_time_nonuniform(y, x):
	pass


def optimise(method="grad_descent", y0=None, plot_time=False,
			 nit = 100, base_lr=1e-3, delta = 1e-3):
	"""Run gradient descent optimisation to optimise Y towards the brachistochrone.
	y0 = initial guess. if none: do linear"""
	Y = np.zeros((nit + 1, X.size))  # the ith element of y gives that iteration's best version of a brachistochrone

	if y0 is None:
		Y[0] = X * H / L  # initialise as a straight line between (0, H) and (L,0)
	else:
		Y[0] = y0

	E = np.eye(X.size) * delta
	deltas = np.ones_like(X) * delta

	lr = base_lr

	f = lambda y: compute_time_uniform(y, X)
	get_grad = lambda y, e: (f(y + eps) - f(y)) / e.max()

	iterations = np.arange(0, nit)
	ts = [] # store performance
	with tqdm(iterations) as tqdm_it:
		for n in tqdm_it:
			y = Y[n] # initial y(x)

			grad = np.zeros_like(X)
			# calc grad for every node except the last (keep this fixed)
			for i in range(1, X.size - 1):
				eps = E[i]
				grad[i] = get_grad(y, eps)

			if "grad_descent" in method:
				change = - lr * grad

			elif method == "newton":
				# Hessian is sparse, with only diagonals and left/right of diagonals being non zero
				hess = np.zeros((X.size, X.size))
				for i in range(X.size):
					eps = E[i]
					hess[i, i] = ( get_grad(y + eps, eps) - get_grad(y, eps) ) / eps.max()
					if i > 0:
						j = i-1
						# finite difference for cross terms
						numerator = f(y + E[i] + E[j]) - f(y + E[j]) - f(y + E[i]) + f(y)
						hess[i, j] = hess[j, i] = numerator / (2 * deltas[i] * deltas[j])

				change = - np.matmul(np.linalg.inv(hess), grad)

			change[[0, -1]] = 0 # no changes allowed for endpoints as they are fixed
			y_new = y + change
			y_new[y_new < 0] = 1e-5  # clamp all values at just above non-zero, to maintain numerical stability

			Y[n+1] = y_new

			t = compute_time_uniform(Y[n + 1], X)
			ts.append(t)
			tqdm_it.set_description(f"Time = {t:.3f}")

	if plot_time:
		fig, tax = plt.subplots()
		tax.plot(ts)

	return Y


def brach(x, t):
	"""Functional to be optimised. t indicates the process of training"""

	if isinstance(x, float) or isinstance(x, int):
		x_idx = np.argmin(np.abs(X - x))

	elif isinstance(x, np.ndarray):
		x_idx = np.array([np.argmin(np.abs(X - i)) for i in x])

	t_idx = np.argmin(np.abs(T-t))

	return Y[t_idx, x_idx]

if __name__ == "__main__":

	method = "newton"
	# vis params
	anim_time = 1  # animation time in seconds
	nit = hyperparams[method]["nit"]
	T = np.linspace(0, anim_time, nit)

	fig, ax = plt.subplots(FigureClass=FuncFig, animate=True, time = anim_time)

	# plot analytical solution
	x_par, y_par = get_brach_x_y(L, H)
	y_res = linear_resample(y_par, x_par, X)
	ax.plot(X, y_res, color="green", lw = 15, alpha = 0.5)

	# plot analytical solution - mu = .5
	x_par, y_par = get_brach_x_y(L, H, mu=0.5)
	y_res = linear_resample(y_par, x_par, X)
	ax.plot(X, y_res, color="red", lw = 15, alpha = 0.5)

	# plt.show()

	# print analytical time results
	brach_analytical_soln = get_brach_times(L, H, X)
	print(f"brach, analytical: T = {brach_analytical_soln[-1]:.4f}")
	print("brach: T = ", compute_time_uniform(y_res, X))

	# Get optimised solution
	Y = optimise(**hyperparams[method])
	print("optimised: T = ", compute_time_uniform(Y[-1], X))


	f = ax.add_func(function=brach,
					t0=T.min(), t1=T.max(), t_samples=fig.fps,
					xrange = X,
					is_functional=True, marker="s")

	pad = 0.1
	strings = [f"Iteration\n{int(n)}" for n in np.linspace(1, nit, fig.n_frames)]
	t = ax.add_text("TEXT", (1-pad)*L, pad*H, strings=strings,
					horizontalalignment="center", fontsize=20,
					bbox=dict(facecolor="blue", alpha=0.5, boxstyle="round"),
					fontfamily=font)  # produce text object


	fig.render()

	# ax formatting
	ax.set_ylim((1+pad)*H, -pad*H)
	ax.set_xlim(-pad*L, (1+pad)*L)
	ax.axis("off")

	fig.subplots_adjust(left=0,bottom=0,top=1,right=1)

	# fig.save(title=f"brach_{method}", fmt="gif")
	plt.show()
