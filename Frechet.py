import numpy as np
import matplotlib.pyplot as plt

def frechet_pdf(z, T, theta):
    z = np.asarray(z)
    pdf = T * theta * z**(-theta-1) * np.exp(-T * z**(-theta))
    pdf[z <= 0] = 0
    return pdf

def frechet_cdf(z, T, theta):
    z = np.asarray(z)
    cdf = np.exp(-T * z**(-theta))
    cdf[z <= 0] = 0
    return cdf

def plot_frechet_distribution(T_values, theta_values, 
                               z_range=(0.1, 5), num_points=1000):   
    z = np.linspace(z_range[0], z_range[1], num_points)
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    fig.suptitle('Frechet Distribution: Productivity in EK Model', fontsize=14, fontweight='bold')
    
    ax1 = axes[0]
    fixed_theta = 4
    for T in T_values:
        pdf = frechet_pdf(z, T, fixed_theta)
        ax1.plot(z, pdf, linewidth=2, label=f'T = {T}')
    ax1.set_title(f'Fixed θ = {fixed_theta}, Varying T', fontsize=11)
    ax1.set_xlabel('Efficiency Level z', fontsize=10)
    ax1.set_ylabel('Probability Density', fontsize=10)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(z_range)
    
    ax2 = axes[1]
    fixed_T = 2
    for theta in theta_values:
        pdf = frechet_pdf(z, fixed_T, theta)
        ax2.plot(z, pdf, linewidth=2, label=f'θ = {theta}')
    ax2.set_title(f'Fixed T = {fixed_T}, Varying θ', fontsize=11)
    ax2.set_xlabel('Efficiency Level z', fontsize=10)
    ax2.set_ylabel('Probability Density', fontsize=10)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(z_range)
    
    plt.tight_layout()
    plt.show()

def interactive_plot():
    # Interactive plotting function (run in Jupyter environment)
    from ipywidgets import interact, FloatSlider
    
    def update_plot(T=2.0, theta=4.0):
        z = np.linspace(0.1, 5, 1000)
        pdf = frechet_pdf(z, T, theta)
        cdf = frechet_cdf(z, T, theta)
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
        
        ax1.plot(z, pdf, 'b-', linewidth=2)
        ax1.set_title(f'Frechet PDF (T={T}, θ={theta})')
        ax1.set_xlabel('Efficiency z')
        ax1.set_ylabel('Probability Density')
        ax1.grid(True, alpha=0.3)
        
        ax2.plot(z, cdf, 'r-', linewidth=2)
        ax2.set_title(f'Frechet CDF (T={T}, θ={theta})')
        ax2.set_xlabel('Efficiency z')
        ax2.set_ylabel('Cumulative Probability')
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()
    
    interact(update_plot,
                T=FloatSlider(min=0.5, max=10.0, step=0.5, value=2.0,
                            description='T:'),
                theta=FloatSlider(min=1, max=20.0, step=0.5, value=4.0,
                                description='θ:'))

if __name__ == "__main__":
    # Example parameters
    T_example = [0.5, 1, 2, 4]
    theta_example = [2, 4, 6, 10]
    
    plot_frechet_distribution(T_values=T_example, theta_values=theta_example)