#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 10:53:51 2020

@author: cfse2
"""
import numpy as np
from framat import Model                  # Import the FramAT model object

model = Model.from_example('cantilever')  # Load the cantilever example model
results = model.run()                               # Run the pre-defined analysis
np.set_printoptions(precision=3)
print(results.get('tensors').get('K'))
print("\n")
print(results.get('tensors').get('M'))
print("\n")
print(np.min(results.get('tensors').get('U')))