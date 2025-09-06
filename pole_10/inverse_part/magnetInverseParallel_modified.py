import deepxde as dde
import numpy as np
import torch
import matplotlib.pyplot as plt
import h5py

d = h5py.File(
    './newProjects/pinn_cylindrical_single_plane_p10_modify_pro/data_generation/H_real.mat')
H_1_observe = ((np.array(d['H_real'])).astype(np.float32)).T

# dde.config.set_default_float("float64")
device = torch.device("cuda:0")

"""defining functions"""


def generate_points(num_rho, num_theta, rho_start, rho_out):
    # This program generates the discrete points INSIDE the magnet!! Not the observation!!!
    # in 2D cases, volume reduces to surface, surface reduces to line
    total_surface_element = num_rho * num_theta
    total_line_element = 2 * num_theta
    rho_pos_ele, rho_neg_ele = (np.zeros((total_surface_element + total_line_element, 1)),
                                np.zeros((total_surface_element + total_line_element, 1)))
    surface_ele = np.zeros((total_surface_element + total_line_element, 1))  # masks

    d_rho = (rho_out - rho_start) / num_rho
    d_theta = 2 * np.pi / num_theta

    totalElement = total_surface_element + total_line_element
    mag_coordinate = np.zeros((totalElement, 2))

    # let's make it simple: fill the first half of 'mag_coordinate' with positive magnetization
    for idx_p in range(1, p + 1, 1):
        count = 0
        starter_positive = (idx_p - 1) * totalElement / (2 * p)
        for ind_theta in range(int((idx_p - 1) * n_theta / p), int((idx_p - 1) * n_theta / p + n_theta / (2 * p)), 1):
            for ind_rho in range(0, n_rho + 2, 1):
                if ind_rho == 0:
                    # line points: inner boundary
                    mag_coordinate[count + int(starter_positive), 0] = rho_start
                    mag_coordinate[count + int(starter_positive), 1] = 0.5 * d_theta + ind_theta * d_theta
                    rho_neg_ele[count + int(starter_positive)] = 1
                    count = count + 1
                elif ind_rho == n_rho + 1:
                    # line points: outer boundary
                    mag_coordinate[count + int(starter_positive), 0] = rho_out
                    mag_coordinate[count + int(starter_positive), 1] = 0.5 * d_theta + ind_theta * d_theta
                    rho_pos_ele[count + int(starter_positive)] = 1
                    count = count + 1
                else:
                    # surface points
                    mag_coordinate[count + int(starter_positive), 0] = rho_start + 0.5 * d_rho + (ind_rho - 1) * d_rho
                    mag_coordinate[count + int(starter_positive), 1] = 0.5 * d_theta + ind_theta * d_theta
                    surface_ele[count + int(starter_positive)] = 1
                    count = count + 1

    # let's make it simple: fill the next half of 'mag_coordinate' with negative magnetization
    for idx_p in range(1, p + 1, 1):
        count = 0
        starter_negative = totalElement / 2 + (idx_p - 1) * totalElement / (2 * p)
        for ind_theta in range(int((idx_p - 1) * n_theta / p + n_theta / (2 * p)),
                               int((idx_p - 1) * n_theta / p + n_theta / p), 1):
            for ind_rho in range(0, n_rho + 2, 1):
                if ind_rho == 0:
                    # line points: inner boundary
                    mag_coordinate[count + int(starter_negative), 0] = rho_start
                    mag_coordinate[count + int(starter_negative), 1] = 0.5 * d_theta + ind_theta * d_theta
                    rho_neg_ele[count + int(starter_negative)] = 1
                    count = count + 1
                elif ind_rho == n_rho + 1:
                    # line points: outer boundary
                    mag_coordinate[count + int(starter_negative), 0] = rho_out
                    mag_coordinate[count + int(starter_negative), 1] = 0.5 * d_theta + ind_theta * d_theta
                    rho_pos_ele[count + int(starter_negative)] = 1
                    count = count + 1
                else:
                    # surface points
                    mag_coordinate[count + int(starter_negative), 0] = rho_start + 0.5 * d_rho + (ind_rho - 1) * d_rho
                    mag_coordinate[count + int(starter_negative), 1] = 0.5 * d_theta + ind_theta * d_theta
                    surface_ele[count + int(starter_negative)] = 1
                    count = count + 1

    # # now, loop all elements (inner+boundary)
    # count = 0
    # for ind_theta in range(0, n_theta, 1):
    #     for ind_rho in range(0, n_rho+2, 1):
    #         if ind_rho == 0:
    #             # line points: inner boundary
    #             mag_coordinate[count, 0] = rho_start
    #             mag_coordinate[count, 1] = 0.5 * d_theta + ind_theta * d_theta
    #             rho_neg_ele[count] = 1
    #             count = count + 1
    #         elif ind_rho == n_rho+1:
    #             # line points: outer boundary
    #             mag_coordinate[count, 0] = rho_out
    #             mag_coordinate[count, 1] = 0.5 * d_theta + ind_theta * d_theta
    #             rho_pos_ele[count] = 1
    #             count = count + 1

    rho_pos_ele = torch.from_numpy(rho_pos_ele).to(device)
    rho_neg_ele = torch.from_numpy(rho_neg_ele).to(device)
    surface_ele = torch.from_numpy(surface_ele).to(device)

    dsurface = d_rho * d_theta * mag_coordinate[:, 0:1]
    dline = d_theta * mag_coordinate[:, 0:1]

    dsurface = torch.from_numpy(dsurface).to(device)
    dline = torch.from_numpy(dline).to(device)

    return mag_coordinate, rho_pos_ele, \
        rho_neg_ele, surface_ele, dsurface, dline


def generate_magnetization(magnet_discrete_coordinate):
    # input: the coordinates of discretization points
    # output: normalized M on the discretization points

    M = np.zeros_like(magnet_discrete_coordinate)
    totalElement = totalSurfaceElement + totalLineElement
    mag = 1

    for idx in range(totalElement):
        if idx < totalElement / 2:
            M[idx, 0] = mag
            M[idx, 1] = 0
        else:
            M[idx, 0] = -mag
            M[idx, 1] = 0

    return M


def func(x):
    # input: the coordinates of discretization points; sigma; x0; y0; z0
    # ATTENTION: do not forget to treat the boundary points at the end the input x.
    M = np.zeros_like(x)
    totalElement = totalSurfaceElement + totalLineElement
    mag = 1

    for idx in range(totalElement):
        if idx < totalElement / 2:
            M[idx, 0] = mag
            M[idx, 1] = 0
        else:
            M[idx, 0] = -mag
            M[idx, 1] = 0

    return M


def pde(x, y):
    ### Hidden codes, contact the author to unlock ###

    return [residual_analytical, residual_M]


def pde_predict(x, y):
    ### Hidden codes, contact the author to unlock ###

    return observeH_sim


"""start main here"""
# scaling coefficients
L_star = 0.01
M_star = 1e6
phi_star = L_star * M_star

# magnet geometry in real dimension
rhoOutReal, rhoInnerReal, thetaLengthReal = 0.01, 0.005, 2 * np.pi
# magnet scaled geometry
rhoOut, rhoInner, thetaLength = rhoOutReal / L_star, rhoInnerReal / L_star, thetaLengthReal
p = 10  # number of pole pairs

# define geometry and sample points of the magnet
n_rho, n_theta = 10, 320
totalSurfaceElement = n_rho * n_theta
totalLineElement = n_theta * 2

# here, the out_element and inner_element are "masks"
magnet_coordinate, out_element, \
    inner_element, surface_element, d_surface, d_line = \
    generate_points(n_rho, n_theta, rhoInner, rhoOut)

geom = dde.geometry.geometry_2d.Disk([0, 0], 1)  # dummy geometry
data = dde.data.PDE(geom, pde, [], num_domain=0, num_boundary=0, anchors=magnet_coordinate, solution=func)

# region   Here, we should import the observation data !!
# Surface 1
H_1_rho_field, H_1_theta_field = H_1_observe[:, :, 0].T, H_1_observe[:, :, 1].T
H_1_rho_reshape, H_1_theta_reshape = (np.reshape(H_1_rho_field, (-1, 1)), np.reshape(H_1_theta_field, (-1, 1)))

H_1_rho_reshape_n = H_1_rho_reshape
H_1_theta_reshape_n = H_1_theta_reshape
H_1_observe = np.concatenate([H_1_rho_reshape_n, H_1_theta_reshape_n], axis=1)
H_1_observe = torch.from_numpy(H_1_observe).to(device)

H_observe = H_1_observe
# endregion

# region Pre-calculation: All observation locations for phi_m, including observeXDifCoord and observeYDifCoord
num_grid = 2000
height_real = 0.002  # in meter

height = height_real / L_star  # note that this is normalized value
delta_h_theta = (2 * np.pi) / num_grid  # note that this is normalized value
delta_h_rho = 0.01  # artificial

observeThetaDifCoord_1_pre = torch.zeros(((num_grid + 1) * 1, 2), device=torch.device("cuda:0"))
observeRhoDifCoord_1_Pos_pre = torch.zeros((num_grid * 1, 2), device=torch.device("cuda:0"))
observeRhoDifCoord_1_Neg_pre = torch.zeros((num_grid * 1, 2), device=torch.device("cuda:0"))

# Processing the rho component
for id_rho in range(0, num_grid, 1):
    observeRhoDifCoord_1_Pos_pre[id_rho, 0] = rhoOut + height + delta_h_rho / 2
    observeRhoDifCoord_1_Neg_pre[id_rho, 0] = rhoOut + height - delta_h_rho / 2

    observeRhoDifCoord_1_Pos_pre[id_rho, 1] = id_rho * delta_h_theta + 0.5 * delta_h_theta
    observeRhoDifCoord_1_Neg_pre[id_rho, 1] = id_rho * delta_h_theta + 0.5 * delta_h_theta

# Processing the theta component
for id_theta in range(0, (num_grid + 1), 1):
    observeThetaDifCoord_1_pre[id_theta, 0] = rhoOut + height
    observeThetaDifCoord_1_pre[id_theta, 1] = id_theta * delta_h_theta
# endregion

# region Here, we generate the distance matrix
# All observation locations (grids) for phi_m
magnet_coordinate = torch.from_numpy(magnet_coordinate).to(device)
# Calculating the distance matrices of each grid

# Processing of observeRhoDifCoord_1_Pos_pre
obs_rhoRhoPos = torch.transpose(torch.unsqueeze(observeRhoDifCoord_1_Pos_pre[:, 0], 1), 0, 1)
obs_rho_RhoPos = obs_rhoRhoPos.repeat(magnet_coordinate[:, 0].size(0), 1)
ele_rhoRhoPos = torch.unsqueeze(magnet_coordinate[:, 0], 1)
ele_rho_RhoPos = ele_rhoRhoPos.repeat(1, observeRhoDifCoord_1_Pos_pre[:, 0].size(0))
obs_thetaRhoPos = torch.transpose(torch.unsqueeze(observeRhoDifCoord_1_Pos_pre[:, 1], 1), 0, 1)
obs_theta_RhoPos = obs_thetaRhoPos.repeat(magnet_coordinate[:, 1].size(0), 1)
ele_thetaRhoPos = torch.unsqueeze(magnet_coordinate[:, 1], 1)
ele_theta_RhoPos = ele_thetaRhoPos.repeat(1, observeRhoDifCoord_1_Pos_pre[:, 1].size(0))
# distance matrix of observeRhoDifCoord_1_Pos_pre
observeRhoDifCoordPos_1_dis = torch.pow((torch.pow(obs_rho_RhoPos, 2) + torch.pow(ele_rho_RhoPos, 2)
                                         - 2 * obs_rho_RhoPos * ele_rho_RhoPos * torch.cos(
            obs_theta_RhoPos - ele_theta_RhoPos)), 0.5)

# Processing of observeRhoDifCoord_1_Neg_pre
obs_rhoRhoNeg = torch.transpose(torch.unsqueeze(observeRhoDifCoord_1_Neg_pre[:, 0], 1), 0, 1)
obs_rho_RhoNeg = obs_rhoRhoNeg.repeat(magnet_coordinate[:, 0].size(0), 1)
ele_rhoRhoNeg = torch.unsqueeze(magnet_coordinate[:, 0], 1)
ele_rho_RhoNeg = ele_rhoRhoNeg.repeat(1, observeRhoDifCoord_1_Neg_pre[:, 0].size(0))

obs_thetaRhoNeg = torch.transpose(torch.unsqueeze(observeRhoDifCoord_1_Neg_pre[:, 1], 1), 0, 1)
obs_theta_RhoNeg = obs_thetaRhoNeg.repeat(magnet_coordinate[:, 1].size(0), 1)
ele_thetaRhoNeg = torch.unsqueeze(magnet_coordinate[:, 1], 1)
ele_theta_RhoNeg = ele_thetaRhoNeg.repeat(1, observeRhoDifCoord_1_Neg_pre[:, 1].size(0))
# distance matrix of observeRhoDifCoord_1_Neg_pre
observeRhoDifCoordNeg_1_dis = torch.pow((torch.pow(obs_rho_RhoNeg, 2) + torch.pow(ele_rho_RhoNeg, 2)
                                         - 2 * obs_rho_RhoNeg * ele_rho_RhoNeg * torch.cos(
            obs_theta_RhoNeg - ele_theta_RhoNeg)), 0.5)

# Processing of "observeThetaDifCoord_1"
obs_rhoTheta = torch.transpose(torch.unsqueeze(observeThetaDifCoord_1_pre[:, 0], 1), 0, 1)
obs_rho_Theta = obs_rhoTheta.repeat(magnet_coordinate[:, 0].size(0), 1)
ele_rhoTheta = torch.unsqueeze(magnet_coordinate[:, 0], 1)
ele_rho_Theta = ele_rhoTheta.repeat(1, observeThetaDifCoord_1_pre[:, 0].size(0))

obs_thetaTheta = torch.transpose(torch.unsqueeze(observeThetaDifCoord_1_pre[:, 1], 1), 0, 1)
obs_theta_Theta = obs_thetaTheta.repeat(magnet_coordinate[:, 1].size(0), 1)
ele_thetaTheta = torch.unsqueeze(magnet_coordinate[:, 1], 1)
ele_theta_Theta = ele_thetaTheta.repeat(1, observeThetaDifCoord_1_pre[:, 1].size(0))

# distance matrix of observeThetaDifCoord_1_pre
observeThetaDifCoord_1_dis = torch.pow((torch.pow(obs_rho_Theta, 2) + torch.pow(ele_rho_Theta, 2)
                                        - 2 * obs_rho_Theta * ele_rho_Theta * torch.cos(
            obs_theta_Theta - ele_theta_Theta)), 0.5)
# endregion

# define network
# layer_size = [2] + [[200]*2, [400]*2, [200]*2] + [2]
layer_size = [2] + [300] + [2]  # If we use FNN_resnet, then these values are "layer_sizes_key", not real layer sizes
activation = 'tanh'  # choose from: ELU, GELU, ReLU, SELU, Sigmoid, SiLU, sin, Swish, tanh
initializer = 'Glorot normal'
num_resnet = 5
net = dde.nn.FNN_resnet(layer_size, num_resnet, activation, initializer)
model = dde.Model(data, net)

# start training
save_model = dde.callbacks.ModelCheckpoint(
    "./newProjects/pinn_cylindrical_single_plane_p10_modify_pro/inverse_part/solve_M.ckpt", period=100000)
model.compile("adam",  # L-BFGS, L-BFGS-B, sgd, rmsprop, adam, adamw
              lr=2e-4,
              loss='mse',
              decay=("step", 2300, 0.75),
              metrics=["l2 relative error"],
              loss_weights=[1, 1])
losshistory, train_state = model.train(iterations=50000,
                                       batch_size=None,
                                       display_every=100,
                                       disregard_previous_best=False,
                                       callbacks=[save_model],
                                       model_restore_path=None,
                                       model_save_path="./newProjects/pinn_cylindrical_single_plane_p10_modify_pro/inverse_part/solve_M.ckpt")

# dde.saveplot(losshistory, train_state, issave=True, isplot=True, output_dir="/home/eegongzhi/PycharmProjects/deepxde/"
#                                                                             "MyProjects/nonuniform_0905/test1_003_noisy/SNR_20/")

# # otherwise, we load the stored model. Disable one of "training" and "restore".
# model.restore(r"./newProjects/pinn_cylindrical_single_plane_p10_modify_pro/inverse_part/solve_M.ckpt-50000.pt", verbose=1)


"""Here, we perform inference"""
magnet_coordinate = magnet_coordinate.detach().cpu().numpy()
M_pred = model.predict(magnet_coordinate) * M_star  # This is the real M
H_pred = model.predict(magnet_coordinate, operator=pde_predict)  # This is the real H
np.savetxt("./newProjects/pinn_cylindrical_single_plane_p10_modify_pro/inverse_part/predicted_M.txt", M_pred)
np.savetxt('./newProjects/pinn_cylindrical_single_plane_p10_modify_pro/inverse_part/predicted_H.txt', H_pred)

# region Error calculation and plotting
# calculate the error of H
H_observe = H_observe.cpu().numpy()
# difference = np.linalg.norm(H_observe - H_pred) / np.linalg.norm(H_observe)
difference = dde.metrics.l2_relative_error(H_observe, H_pred)
print("Relative error of H = ", difference)

# calculate the error of M: only for simple M. For complicated M, we should import this from separate files!!!!
M_true = M_star * generate_magnetization(magnet_coordinate)  # all in scaled value
difference = dde.metrics.l2_relative_error(M_true, M_pred)
print("Relative error of M = ", difference)

"""Plot the observation H and the predicted H"""
observeLocationTheta = np.zeros((num_grid * 1, 1))
for id_theta in range(0, num_grid, 1):
    observeLocationTheta[id_theta, 0] = id_theta * delta_h_theta + 0.5 * delta_h_theta

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(observeLocationTheta, H_observe[:, 0], '-', label='Observed H_rho')
ax.plot(observeLocationTheta, H_pred[:, 0], '--', label='Simulated H_rho')
ax2 = ax.twinx()
ax.plot(observeLocationTheta, H_observe[:, 1], '-', label='Observed H_theta')
ax.plot(observeLocationTheta, H_pred[:, 1], '--', label='Simulated H_theta')

ax.legend(loc=0)
ax.grid()
ax.set_xlabel("Theta (rad)")
ax.set_ylabel(r"H_rho")
ax2.set_ylabel(r"H_theta")

plt.show()

