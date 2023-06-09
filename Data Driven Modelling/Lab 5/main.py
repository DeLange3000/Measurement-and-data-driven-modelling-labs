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
import numpy as num


dataset = pd.read_csv('household_power_consumption.csv')
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

# Select additional columns as input features
additional_columns = ['Global_reactive_power','Voltage','Global_intensity']
X = dataset[additional_columns] # input data set
y = dataset['Global_active_power']  # Target variable (ouput data set)
X = torch.tensor(X.values,dtype=torch.float32)
y = torch.tensor(y.values, dtype=torch.float32)

# First split and then normalize (If we do the inverse, in real word it doesn't match, because normalise data took in to account the mean and the variance, and it will make worng the data during the standardization)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42) #test size the percent of data that we have to preprocess

# Normalize the input features & output
scaler = MinMaxScaler()
X_train = scaler.fit_transform(X_train) #normalise data
X_train = torch.tensor(X_train,dtype=torch.float32)
X_test = scaler.fit_transform(X_test)
X_test = torch.tensor(X_test,dtype=torch.float32)
y_train = y_train.view(y_train.shape[0], 1) #change output tensor layout
y_test = y_test.view(y_test.shape[0], 1)

class NeuralNetwork(nn.Module):
    def __init__(self):
        super(NeuralNetwork, self).__init__()
        self.flatten = nn.Flatten()
        self.layers = nn.Sequential(
            nn.Linear(input_size, hidden_size), #input layer -> hidden layer
            nn.ReLU(), #activation function
            # nn.Linear(hidden_size, hidden_size),
            # nn.ReLU(),
            # nn.Linear(hidden_size, hidden_size),
            # nn.ReLU(),
            nn.Linear(hidden_size, output_size), #hidden layer -> output layer
        )

    def forward(self, x):
        x = self.flatten(x)
        y = self.layers(x)
        return y

input_size = X.size(dim = 1) #set input layer size (amount of input layer variables)
hidden_sizes = [50, 150, 250, 350, 500, 1000] #amount of neurons per layer
output_size = 1

epochs = 300 # Number of iterations of the training loop
for hidden_size in hidden_sizes:

    net = NeuralNetwork().to(device) #create neural network
    print(net)

    # Define loss function and optimizer
    criterion = nn.MSELoss()
    optimizer = optim.SGD(net.parameters(), lr=0.01) #low rate lr (to recheck)

    # Training loop
    test_loss_array = []
    for epoch in range(epochs):
        ### Training
        # Put the model in training mode
        #here is the backward implementation

        net.train()

        y_pred = net(X_train) #output prediction
        loss = criterion(y_pred, y_train) #MSE
        optimizer.zero_grad() # Avoid accumulation of gradients from previous iterations
        loss.backward() # Compute the gradient
        optimizer.step() # Set the new weights of the model as a function of the gradient

        ### Testing
        # Put the model in evaluation mode
        net.eval()

        with torch.inference_mode():
            test_pred = net(X_test)
            test_loss = criterion(test_pred, y_test)

        test_loss_array.append(test_loss) #get MSE from testing
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

    x = num.arange(len(test_loss_array))

    line = plt.plot(x, test_loss_array, label = str(hidden_size))
plt.xlabel('Epochs')
plt.ylabel('Mean Square Error')
plt.legend(hidden_sizes)
plt.show()
