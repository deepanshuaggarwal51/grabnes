import numpy as np
import matplotlib.pyplot as plt

def parse_bz_file(filename):
    """
    Parse a single Brillouin zone data file
    """
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Warning: {filename} not found, skipping this layer")
        return None, None, None, None
    
    # Parse the data
    bz_vertices = []
    k_point = None
    b1 = None
    b2 = None
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        # Skip comments and empty lines
        if line.startswith('#') or not line:
            # Check for specific section headers
            if '# BZ vertices (x, y):' in line:
                # Read next 6 lines for BZ vertices
                for j in range(1, 7):
                    if i + j < len(lines):
                        coords = list(map(float, lines[i + j].strip().split()))
                        bz_vertices.append(coords)
                i += 6
            elif '# K-point:' in line:
                # Read next line for K-point
                if i + 1 < len(lines):
                    k_point = list(map(float, lines[i + 1].strip().split()))
                i += 1
            elif '# Reciprocal lattice vectors b1, b2:' in line:
                # Read next 2 lines for b1 and b2
                if i + 1 < len(lines):
                    b1 = list(map(float, lines[i + 1].strip().split()))
                if i + 2 < len(lines):
                    b2 = list(map(float, lines[i + 2].strip().split()))
                i += 2
        i += 1
    
    # Convert to numpy arrays
    bz_vertices = np.array(bz_vertices) if bz_vertices else None
    k_point = np.array(k_point) if k_point else None
    b1 = np.array(b1) if b1 else None
    b2 = np.array(b2) if b2 else None
    
    return bz_vertices, k_point, b1, b2

def plot_three_layer_brillouin_zones(filenames=None):
    """
    Plot the Brillouin zones for three layers from separate files
    """
    if filenames is None:
        filenames = [
            'brillouin_zones_debug_K1.dat',
            'brillouin_zones_debug_K2.dat', 
            'brillouin_zones_debug_K3.dat'
        ]
    
    # Colors and labels for the three layers
    colors = ['blue', 'red', 'green']
    fill_colors = ['lightblue', 'lightcoral', 'lightgreen']
    labels = ['Layer 1', 'Layer 2', 'Layer 3']
    k_markers = ['*', 's', '^']  # star, square, triangle
    
    # Create the plot
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    
    # Parse and plot each layer
    for i, (filename, color, fill_color, label, k_marker) in enumerate(zip(filenames, colors, fill_colors, labels, k_markers)):
        bz_vertices, k_point, b1, b2 = parse_bz_file(filename)
        
        if bz_vertices is None:
            continue
            
        # Plot the Brillouin zone (close the hexagon by adding first point at end)
        if len(bz_vertices) == 6:
            bz_closed = np.vstack([bz_vertices, bz_vertices[0]])
            ax.plot(bz_closed[:, 0], bz_closed[:, 1], color=color, linewidth=2, label=f'{label} BZ')
            ax.fill(bz_closed[:, 0], bz_closed[:, 1], alpha=0.15, color=fill_color)
            
            # Plot BZ vertices
            ax.scatter(bz_vertices[:, 0], bz_vertices[:, 1], c=color, s=30, zorder=5, alpha=0.7)
        
        # Plot K-point
        if k_point is not None:
            ax.scatter(k_point[0], k_point[1], c=color, s=120, marker=k_marker, 
                      zorder=6, label=f'{label} K-point', edgecolors='black', linewidth=1)
        
        # Plot reciprocal lattice vectors from origin (only for first layer to avoid clutter)
        if i == 0 and b1 is not None and b2 is not None:
            arrow_scale = 0.8  # Make arrows a bit smaller
            head_size = 0.03
            ax.arrow(0, 0, b1[0]*arrow_scale, b1[1]*arrow_scale, 
                    head_width=head_size, head_length=head_size, fc='darkgreen', ec='darkgreen', 
                    alpha=0.7, zorder=3, label='b1')
            ax.arrow(0, 0, b2[0]*arrow_scale, b2[1]*arrow_scale, 
                    head_width=head_size, head_length=head_size, fc='darkorange', ec='darkorange', 
                    alpha=0.7, zorder=3, label='b2')
            
            # Add text labels for vectors
            ax.text(b1[0]*arrow_scale*1.1, b1[1]*arrow_scale*1.1, 'b1', fontsize=11, 
                   color='darkgreen', fontweight='bold')
            ax.text(b2[0]*arrow_scale*1.1, b2[1]*arrow_scale*1.1, 'b2', fontsize=11, 
                   color='darkorange', fontweight='bold')
        
        # Print layer info
        if k_point is not None:
            print(f"{label}:")
            print(f"  K-point: ({k_point[0]:.6f}, {k_point[1]:.6f})")
            print(f"  |K| = {np.linalg.norm(k_point):.6f}")
            
        if b1 is not None and b2 is not None:
            print(f"  |b1| = {np.linalg.norm(b1):.6f}")
            print(f"  |b2| = {np.linalg.norm(b2):.6f}")
            
            # Check if vectors are 60 degrees apart (should be for graphene)
            dot_product = np.dot(b1, b2)
            magnitudes = np.linalg.norm(b1) * np.linalg.norm(b2)
            if magnitudes > 0:
                angle = np.arccos(np.clip(dot_product / magnitudes, -1, 1)) * 180 / np.pi
                print(f"  Angle between b1 and b2: {angle:.1f}°")
        print()
    
    # Add origin
    ax.scatter(0, 0, c='black', s=80, marker='o', zorder=5, label='Origin')
    
    # Load and plot G-vectors with error handling
    try:
        gvecs = np.loadtxt('g_vectors_debug.dat')
        print(f"Loaded {len(gvecs)} G-vectors")
        print(f"G-vector range: x=[{gvecs[:, 0].min():.3f}, {gvecs[:, 0].max():.3f}], y=[{gvecs[:, 1].min():.3f}, {gvecs[:, 1].max():.3f}]")
        
        # Check if G-vectors are in reasonable range compared to BZ
        if len(gvecs) > 0:
            g_max = max(abs(gvecs[:, 0].max()), abs(gvecs[:, 0].min()), 
                       abs(gvecs[:, 1].max()), abs(gvecs[:, 1].min()))
            print(f"Max G-vector magnitude: {g_max:.3f}")
            
            # Plot G-vectors with a neutral color that works with all three layers
            ax.scatter(gvecs[:, 0], gvecs[:, 1], c='gray', s=15, alpha=0.5, 
                      label=f'G-vectors ({len(gvecs)})', zorder=2)
            
            # Print first few G-vectors for debugging
            print("First 5 G-vectors:")
            for i in range(min(5, len(gvecs))):
                print(f"  G({i+1}): ({gvecs[i, 0]:.6f}, {gvecs[i, 1]:.6f})")
        else:
            print("No G-vectors found in file")
            
    except FileNotFoundError:
        print("Warning: g_vectors_debug.dat not found")
    except Exception as e:
        print(f"Error loading G-vectors: {e}")
    
    # Set equal aspect ratio and grid
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Labels and title
    ax.set_xlabel('kx (1/Å)', fontsize=12)
    ax.set_ylabel('ky (1/Å)', fontsize=12)
    ax.set_title('Three-Layer Graphene Brillouin Zones', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.show()
    
    return fig, ax

def plot_single_brillouin_zone(filename='brillouin_zones_debug.dat'):
    """
    Plot a single Brillouin zone (backward compatibility)
    """
    bz_vertices, k_point, b1, b2 = parse_bz_file(filename)
    
    if bz_vertices is None:
        print(f"Could not load data from {filename}")
        return None, None
    
    # Create the plot
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    
    # Plot the Brillouin zone (close the hexagon by adding first point at end)
    if len(bz_vertices) == 6:
        bz_closed = np.vstack([bz_vertices, bz_vertices[0]])
        ax.plot(bz_closed[:, 0], bz_closed[:, 1], 'b-', linewidth=2, label='Brillouin Zone')
        ax.fill(bz_closed[:, 0], bz_closed[:, 1], alpha=0.2, color='lightblue')
        
        # Plot BZ vertices
        ax.scatter(bz_vertices[:, 0], bz_vertices[:, 1], c='blue', s=50, zorder=5, label='BZ vertices')
    
    # Plot K-point
    if k_point is not None:
        ax.scatter(k_point[0], k_point[1], c='red', s=100, marker='*', zorder=6, label='K-point')
    
    # Plot reciprocal lattice vectors from origin
    if b1 is not None and b2 is not None:
        ax.arrow(0, 0, b1[0], b1[1], head_width=0.05, head_length=0.05, fc='green', ec='green', label='b1')
        ax.arrow(0, 0, b2[0], b2[1], head_width=0.05, head_length=0.05, fc='orange', ec='orange', label='b2')
        
        # Add text labels for vectors
        ax.text(b1[0]*1.1, b1[1]*1.1, 'b1', fontsize=12, color='green', fontweight='bold')
        ax.text(b2[0]*1.1, b2[1]*1.1, 'b2', fontsize=12, color='orange', fontweight='bold')
    
    # Add origin
    ax.scatter(0, 0, c='black', s=50, marker='o', zorder=5, label='Origin')
    
    # Load and plot G-vectors with error handling
    try:
        gvecs = np.loadtxt('g_vectors_debug.dat')
        print(f"Loaded {len(gvecs)} G-vectors")
        print(f"G-vector range: x=[{gvecs[:, 0].min():.3f}, {gvecs[:, 0].max():.3f}], y=[{gvecs[:, 1].min():.3f}, {gvecs[:, 1].max():.3f}]")
        
        # Check if G-vectors are in reasonable range compared to BZ
        if len(gvecs) > 0:
            g_max = max(abs(gvecs[:, 0].max()), abs(gvecs[:, 0].min()), 
                       abs(gvecs[:, 1].max()), abs(gvecs[:, 1].min()))
            print(f"Max G-vector magnitude: {g_max:.3f}")
            
            # Plot G-vectors
            ax.scatter(gvecs[:, 0], gvecs[:, 1], c='red', s=20, alpha=0.6, 
                      label=f'G-vectors ({len(gvecs)})', zorder=4)
            
            # Print first few G-vectors for debugging
            print("First 5 G-vectors:")
            for i in range(min(5, len(gvecs))):
                print(f"  G({i+1}): ({gvecs[i, 0]:.6f}, {gvecs[i, 1]:.6f})")
        else:
            print("No G-vectors found in file")
            
    except FileNotFoundError:
        print("Warning: g_vectors_debug.dat not found")
    except Exception as e:
        print(f"Error loading G-vectors: {e}")
    
    # Set equal aspect ratio and grid
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # Labels and title
    ax.set_xlabel('kx (1/Å)', fontsize=12)
    ax.set_ylabel('ky (1/Å)', fontsize=12)
    ax.set_title('Graphene Brillouin Zone', fontsize=14, fontweight='bold')
    
    # Print some info
    if k_point is not None:
        print(f"K-point: ({k_point[0]:.6f}, {k_point[1]:.6f})")
        print(f"|K| = {np.linalg.norm(k_point):.6f}")
    
    if b1 is not None and b2 is not None:
        print(f"b1: ({b1[0]:.6f}, {b1[1]:.6f})")
        print(f"b2: ({b2[0]:.6f}, {b2[1]:.6f})")
        print(f"|b1| = {np.linalg.norm(b1):.6f}")
        print(f"|b2| = {np.linalg.norm(b2):.6f}")
        
        # Check if vectors are 60 degrees apart (should be for graphene)
        dot_product = np.dot(b1, b2)
        magnitudes = np.linalg.norm(b1) * np.linalg.norm(b2)
        if magnitudes > 0:
            angle = np.arccos(np.clip(dot_product / magnitudes, -1, 1)) * 180 / np.pi
            print(f"Angle between b1 and b2: {angle:.1f}°")
    
    plt.tight_layout()
    plt.show()
    
    return fig, ax

# Usage examples:
if __name__ == "__main__":
    # Plot three-layer Brillouin zones
    print("Plotting three-layer Brillouin zones...")
    fig, ax = plot_three_layer_brillouin_zones()
    
    # Or plot a single BZ (backward compatibility)
    # fig, ax = plot_single_brillouin_zone()
