# -*- coding: utf-8 -*-
"""DDM-Lab5.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/15fbmtHrqe_My9xjCyt5OLiHnyli-E03o
"""

import torch
from torch import nn
import torch.optim as optim
from torch.utils.data import DataLoader
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import mean_squared_error

from google.colab import drive
drive.mount('/content/drive')
data_path = '/content/drive/MyDrive'

dataset = pd.read_csv(data_path +'/household_power_consumption.csv')
device = "cpu"

"""Source: https://pytorch.org/tutorials/beginner/basics/quickstart_tutorial.html"""

dataset = dataset.drop(columns=['time'])
print(dataset.head())

# # Date and time can also be removed, but they will affect power consumption so it's good to include them
# # Convert the date and time column to a pandas datetime object
# dataset['time'] = pd.to_datetime(dataset['time'])

# # Extract relevant features from the datetime
# dataset['year'] = dataset['time'].dt.year
# dataset['month'] = dataset['time'].dt.month
# dataset['day'] = dataset['time'].dt.day
# dataset['hour'] = dataset['time'].dt.hour
# dataset['minute'] = dataset['time'].dt.minute
# dataset['second'] = dataset['time'].dt.second

dataset = dataset.loc[0:50000]

dataset.shape

# Select additional columns as input features
additional_columns = ['Global_reactive_power','Voltage','Global_intensity']
X = dataset[additional_columns]
y = dataset['Global_active_power']  # Target variable

# Normalize the input features & output
scaler = MinMaxScaler()
X_normalized = scaler.fit_transform(X)

X_tensor = torch.tensor(X_normalized, dtype=torch.float32)
y_tensor = torch.tensor(y.values, dtype=torch.float32)

X_train, X_test, y_train, y_test = train_test_split(X_tensor, y_tensor, test_size=0.2, random_state=42)

class NeuralNetwork(nn.Module):
    def __init__(self):
        super(NeuralNetwork, self).__init__()
        self.layers = nn.Sequential(
          nn.Linear(input_size, hidden_size),
          nn.ReLU(),
          nn.Linear(hidden_size, output_size),
        )

    def forward(self, x):        
        y = self.layers(x)
        return y

input_size = X_tensor.size(dim = 1)
hidden_size = 20
output_size = 1

net = NeuralNetwork().to(device)
print(net)

# Define loss function and optimizer
criterion = nn.MSELoss()
optimizer = optim.Adam(net.parameters(), lr=0.0001)

num_epochs = 10

# Training loop
net.train()

epochs = 10 # Number of iterations of the training loop
for epoch in range(epochs):
    ### Training
    # Put the model in training mode
    net.train()

    y_pred = net(X_train)
    loss = criterion(y_pred, y_train)
    optimizer.zero_grad() # Avoid accumulation of gradients from previous iterations
    loss.backward() # Compute the gradient
    optimizer.step() # Set the new weights of the model as a function of the gradient

    ### Testing
    # Put the model in evaluation mode
    net.eval()

    with torch.inference_mode():
      test_pred = net(X_test)
      test_loss = criterion(test_pred, y_test)

# for epoch in range(num_epochs):
#     ### Training
#     # Put the model in training mode
#     net.train()

#     y_pred = net(X_train)
#     y_train = y_train.type(torch.float32)
#     loss = criterion(y_pred, y_train)
#     optimizer.zero_grad() # Avoid accumulation of gradients from previous iterations
#     loss.backward() # Compute the gradient
#     optimizer.step() # Set the new weights of the model as a function of the gradient

#     ### Testing
#     # Put the model in evaluation mode
#     net.eval()