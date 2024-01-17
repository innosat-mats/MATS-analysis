#%%
import numpy as np
import matplotlib.pyplot as plt

# Define the 1D sinusoidal wave
x = np.linspace(0, 3 * 2 * np.pi, 100)  # Three wavelengths
waveinput = np.sin(x)

# Create the 2D field
#field = np.outer(waveinput, waveinput)

# Create the 2D fields
field_x = np.outer(waveinput, np.ones_like(waveinput))
field_y = np.outer(np.ones_like(waveinput), waveinput)

# Sum the 2D fields
field = field_x + field_y


# Plot the 2D field
plt.imshow(field, cmap='viridis', origin='lower')
plt.colorbar(label='Amplitude')
plt.title('2D Field from 1D Sinusoidal Wave')
plt.show()
# %%
