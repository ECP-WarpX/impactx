#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Ryan Sandberg
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from enum import Enum

import torch
from torch import nn


class Activation(Enum):
    """
    Activation class provides an enumeration type for the supported activation layers
    """

    ReLU = 1
    Tanh = 2
    PReLU = 3
    Sigmoid = 4


def get_enum_type(type_to_test, EnumClass):
    """
    Returns the enumeration type associated to type_to_test in EnumClass

    Parameters
    ----------
    type_to_test: EnumClass, int or str
        object whose Enum class is to be obtained
    EnumClass: Enum class
        Enum class to test
    """
    if type(type_to_test) is EnumClass:  ## Useful ?
        return type_to_test
    if type(type_to_test) is int:
        return EnumClass(type_to_test)
    if type(type_to_test) is str:
        return getattr(EnumClass, type_to_test)


class ConnectedNN(nn.Module):
    """
    ConnectedNN is a class of fully connected neural networks
    """

    def __init__(self, layers):
        super().__init__()
        self.stack = nn.Sequential(*layers)

    def forward(self, x):
        return self.stack(x)


class OneActNN(ConnectedNN):
    """
    OneActNN is class of fully connected neural networks admitting only one activation function
    """

    def __init__(self, n_in, n_out, n_hidden_nodes, n_hidden_layers, act):
        self.n_in = n_in
        self.n_out = n_out
        self.n_hidden_layers = n_hidden_layers
        self.n_hidden_nodes = n_hidden_nodes
        self.act = act

        layers = [nn.Linear(self.n_in, self.n_hidden_nodes)]

        for ii in range(self.n_hidden_layers):
            if self.act is Activation.ReLU:
                layers += [nn.ReLU()]
            if self.act is Activation.Tanh:
                layers += [nn.Tanh()]
            if self.act is Activation.PReLU:
                layers += [nn.PReLU()]
            if self.act is Activation.Sigmoid:
                layers += [nn.Sigmoid()]

            if ii < self.n_hidden_layers - 1:
                layers += [nn.Linear(self.n_hidden_nodes, self.n_hidden_nodes)]

        layers += [nn.Linear(self.n_hidden_nodes, self.n_out)]

        super().__init__(layers)


class surrogate_model:
    """ """

    def __init__(self, dataset_file, model_file):
        self.dataset = torch.load(dataset_file)
        model_dict = torch.load(model_file, map_location=torch.device("cpu"))
        n_in = model_dict["model_state_dict"]["stack.0.weight"].shape[1]
        final_layer_key = list(model_dict["model_state_dict"].keys())[-1]
        n_out = model_dict["model_state_dict"][final_layer_key].shape[0]
        n_hidden_nodes = model_dict["model_state_dict"]["stack.0.weight"].shape[0]
        n_hidden_layers = int(len(model_dict["model_state_dict"].keys()) / 2 - 1)

        activation_type = model_dict["activation"]
        activation = get_enum_type(activation_type, Activation)

        self.neural_network = OneActNN(
            n_in=n_in,
            n_out=n_out,
            n_hidden_nodes=n_hidden_nodes,
            n_hidden_layers=n_hidden_layers,
            act=activation,
        )
        self.neural_network.load_state_dict(model_dict["model_state_dict"])
        self.neural_network.eval()

    def __call__(self, data_arr):
        data_arr -= self.dataset["source_means"]
        data_arr /= self.dataset["source_stds"]
        data_arr = data_arr.float()
        with torch.no_grad():
            data_arr_post_model = self.neural_network(data_arr)

        data_arr_post_model *= self.dataset["target_stds"]
        data_arr_post_model += self.dataset["target_means"]
        return data_arr_post_model
