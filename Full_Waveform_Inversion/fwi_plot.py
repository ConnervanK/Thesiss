import numpy as np
import matplotlib.pyplot as plt

def plot_fwi_results():
    SHAPE = (500, 250)      
    EXTENT = (1.0, 0.5)     

    try:
        m_inverted = np.load('Full_Waveform_Inversion/fwi_result.npy')
        
        # Convert slowness squared back to velocity map
        vel_inverted = 1.0 / np.sqrt(m_inverted)
        
        # Convert to relative permittivity: er = (c / v)^2, with c = 0.3 m/ns
        er_inverted = (0.3 / vel_inverted)**2
        
        plt.figure(figsize=(10, 5))
        plt.imshow(er_inverted.T, extent=[0, EXTENT[0], EXTENT[1], 0],
                   cmap='viridis', vmax=10, vmin=1)
        plt.title('FWI Inverted Relative Permittivity')
        plt.colorbar(label='$\epsilon_r$')
        plt.xlabel('X (m)')
        plt.ylabel('Depth (m)')
        
        print("Saving plot to Full_Waveform_Inversion/fwi_permittivity.png")
        plt.savefig('Full_Waveform_Inversion/fwi_permittivity.png')
        plt.show()
    except FileNotFoundError:
        print("FWI result not found. Run the inversion first.")

if __name__ == '__main__':
    plot_fwi_results()
