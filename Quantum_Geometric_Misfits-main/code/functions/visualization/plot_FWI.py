import matplotlib.pyplot as plt
import numpy as np
import io
import imageio.v2 as imageio

def plot_fwi_results(v_true, v_current, v_init, grad, src_coords, rec_coords, nx, nz, dx, dz, iteration, save_path=None, source_frequency=None):
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f'FWI Iteration {iteration} (Freq {source_frequency} Hz)')
    
    # extent: left, right, bottom, top
    extent = [0, nx*dx, nz*dz, 0]
    
    images = []
    titles = ['True Permittivity ($\\epsilon_r$)', 'Initial Permittivity ($\\epsilon_r$)', 'Current Permittivity ($\\epsilon_r$)', 'Gradient']
    data = [
        v_true.reshape(nx, nz).T, 
        v_init.reshape(nx, nz).T, 
        v_current.reshape(nx, nz).T, 
        grad.reshape(nx, nz).T
    ]
    
    for i, ax in enumerate(axs.ravel()):
        im = ax.imshow(data[i], extent=extent, cmap='jet' if i < 3 else 'RdBu', aspect='auto')
        ax.set_title(titles[i])
        ax.set_xlabel('x (m)')
        ax.set_ylabel('z (m)')
        plt.colorbar(im, ax=ax)
        
        # Plot sources and receivers
        if i < 3:
            if src_coords is not None:
                src_arr = np.array(src_coords)
                if src_arr.ndim == 2:
                    ax.scatter(src_arr[:, 0], src_arr[:, 1], c='red', marker='*', label='Sources', s=50)
            if rec_coords is not None:
                rec_arr = np.array(rec_coords)
                if rec_arr.ndim == 2:
                    ax.scatter(rec_arr[:, 0], rec_arr[:, 1], c='black', marker='v', label='Receivers', s=20)
            if i == 0 and (src_coords is not None or rec_coords is not None):
                ax.legend(loc='upper right')

    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path)
    
    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=100)
    plt.close(fig)
    buf.seek(0)
    img_arr = imageio.imread(buf)
    return img_arr

def plot_convergence(history, save_path=None):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(history['total_loss'], label='Total Loss')
    if 'geometric_misfit' in history and len(history['geometric_misfit']) > 0:
        ax.plot(history['geometric_misfit'], label='Geometric/Data Misfit')
    if 'param_misfit' in history and len(history['param_misfit']) > 0:
        ax.plot(history['param_misfit'], label='Parameter Misfit')
        
    ax.set_yscale('log')
    ax.set_title('FWI Convergence')
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Loss (log scale)')
    ax.legend()
    ax.grid(True)
    
    if save_path:
        plt.savefig(save_path)
    plt.close(fig)