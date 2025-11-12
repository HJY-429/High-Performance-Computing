import matplotlib.pyplot as plt
import numpy as np

def get_norm_num():
    numbers = []
    bias = 3  # since k=3, bias = 2^(3-1) - 1 = 3
    for e in range(1, 7):
        E = e - bias
        for m in range(4):  # mantissa: 00, 01, 10, 11
            mantissa = 1.0 + m * 0.25  # 1.00, 1.25, 1.50, 1.75
            value = mantissa * (2 ** E)
            numbers.append(value)
            numbers.append(-value)
    return numbers

def get_denorm_num():
    numbers = []
    bias = 3
    E = 1 - bias  # -2 for denormalized
    for m in range(4):  # mantissa: 00, 01, 10, 11
        mantissa_value = m * 0.25  # 0.00, 0.25, 0.50, 0.75
        value = mantissa_value * (2 ** E)
        numbers.append(value)
        numbers.append(-value)
    return numbers

def plot_num(normalized, denormalized):
    
    # Plot settings
    plt.figure(figsize=(10, 2))
    plt.title("Representable Numbers in a 6 bit Floating Point System")
    plt.xlabel("Value")
    plt.yticks([])
    
    plt.plot(normalized, np.zeros_like(normalized), 'bo', markersize=5, label='Normalized')
    plt.plot(denormalized, np.zeros_like(denormalized), 'r*', markersize=5, label='Denormalized')
    
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.show()

# Generate numbers
normalized = get_norm_num()
denormalized = get_denorm_num()

# Plot
plot_num(normalized, denormalized)

print("Normalized numbers:", sorted(normalized))
print("Denormalized numbers:", sorted(denormalized))