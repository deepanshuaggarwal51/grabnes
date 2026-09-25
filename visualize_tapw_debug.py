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

def parse_moire_data(filename='moire_bz_debug.dat'):
    """
    Parse moiré BZ data from Fortran output
    """
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Warning: {filename} not found")
        return None, None
    
    moire_vertices = []
    moire_rcell = None
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        if '# Moire BZ vertices (x, y):' in line:
            # Read next 6 lines for BZ vertices
            for j in range(1, 7):
                if i + j < len(lines):
                    coords = list(map(float, lines[i + j].strip().split()))
                    moire_vertices.append(coords)
            i += 6
        elif '# Moire rcell vectors:' in line:
            # Read next 2 lines for rcell vectors
            rcell_data = []
            for j in range(1, 3):
                if i + j < len(lines):
                    coords = list(map(float, lines[i + j].strip().split()[:2]))  # Only x,y
                    rcell_data.append(coords)
            if len(rcell_data) == 2:
                moire_rcell = np.array(rcell_data)
            i += 2
        i += 1
    
    moire_vertices = np.array(moire_vertices) if moire_vertices else None
    
    return moire_vertices, moire_rcell

def parse_kpath_data(filename='kpath_debug.dat'):
    """
    Parse k-path data from Fortran output
    """
    try:
        kpath_frac = np.loadtxt(filename + '_fractional')
        kpath_abs = np.loadtxt(filename + '_absolute')
        return kpath_frac, kpath_abs
    except FileNotFoundError:
        print(f"Warning: k-path files not found")
        return None, None

def plot_comprehensive_tapw_debug():
    """
    Comprehensive TAPW debugging visualization for 2-layer system
    """
    # Create a large figure with multiple subplots
    fig = plt.figure(figsize=(20, 12))
    
    # Define subplot layout: 2x3 grid
    gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)
    
    # Colors for layers
    colors = ['blue', 'red']
    fill_colors = ['lightblue', 'lightcoral']
    labels = ['Layer 1', 'Layer 2']
    
    # ========== SUBPLOT 1: Graphene BZ with G-vectors ==========
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('Graphene BZ + G-vectors', fontweight='bold')
    
    # Load and plot both graphene layers
    filenames = ['brillouin_zones_debug_K1.dat', 'brillouin_zones_debug_K2.dat']
    
    all_k_points = []
    for i, (filename, color, fill_color, label) in enumerate(zip(filenames, colors, fill_colors, labels)):
        bz_vertices, k_point, b1, b2 = parse_bz_file(filename)
        
        if bz_vertices is not None:
            # Plot BZ
            bz_closed = np.vstack([bz_vertices, bz_vertices[0]])
            ax1.plot(bz_closed[:, 0], bz_closed[:, 1], color=color, linewidth=2, label=f'{label} BZ')
            ax1.fill(bz_closed[:, 0], bz_closed[:, 1], alpha=0.15, color=fill_color)
            
            # Store K-point
            if k_point is not None:
                all_k_points.append(k_point)
                ax1.scatter(k_point[0], k_point[1], c=color, s=120, marker='*', 
                           zorder=6, label=f'{label} K-point', edgecolors='black', linewidth=1)
        
        # Plot reciprocal vectors (only first layer)
        if i == 0 and b1 is not None and b2 is not None:
            ax1.arrow(0, 0, b1[0], b1[1], head_width=0.03, head_length=0.03, 
                     fc='darkgreen', ec='darkgreen', alpha=0.7, zorder=3, label='b1')
            ax1.arrow(0, 0, b2[0], b2[1], head_width=0.03, head_length=0.03, 
                     fc='darkorange', ec='darkorange', alpha=0.7, zorder=3, label='b2')
    
    # Add origin and G-vectors
    ax1.scatter(0, 0, c='black', s=80, marker='o', zorder=5, label='Origin')
    
    # Load G-vectors
    try:
        gvecs = np.loadtxt('g_vectors_debug.dat')
        ax1.scatter(gvecs[:, 0], gvecs[:, 1], c='purple', s=15, alpha=0.7, 
                   label=f'G-vectors ({len(gvecs)})', zorder=4)
        print(f"Loaded {len(gvecs)} G-vectors")
    except:
        print("No G-vectors found")
    
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=8)
    ax1.set_xlabel('kx (1/Å)')
    ax1.set_ylabel('ky (1/Å)')
    
    # ========== SUBPLOT 2: Moiré BZ at Origin ==========
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.set_title('Moiré BZ at Origin (0,0)', fontweight='bold')
    
    # Parse moiré BZ data
    moire_vertices, moire_rcell = parse_moire_data()
    
    if moire_vertices is not None:
        # Plot moiré BZ
        moire_closed = np.vstack([moire_vertices, moire_vertices[0]])
        ax2.plot(moire_closed[:, 0], moire_closed[:, 1], 'green', linewidth=2, label='Moiré BZ')
        ax2.fill(moire_closed[:, 0], moire_closed[:, 1], alpha=0.2, color='lightgreen')
        
        # Plot moiré rcell vectors
        if moire_rcell is not None:
            ax2.arrow(0, 0, moire_rcell[0, 0], moire_rcell[0, 1], 
                     head_width=0.002, head_length=0.002, fc='red', ec='red', label='rcell[1]')
            ax2.arrow(0, 0, moire_rcell[1, 0], moire_rcell[1, 1], 
                     head_width=0.002, head_length=0.002, fc='orange', ec='orange', label='rcell[2]')
    
    ax2.scatter(0, 0, c='black', s=80, marker='o', zorder=5, label='Origin')
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=8)
    ax2.set_xlabel('kx (1/Å)')
    ax2.set_ylabel('ky (1/Å)')
    
    # ========== SUBPLOT 3: Moiré BZ around k_ref ==========
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.set_title('Moiré BZ around k_ref', fontweight='bold')
    
    # Plot same moiré BZ but shifted to k_ref location
    if moire_vertices is not None and len(all_k_points) > 0:
        k_ref = all_k_points[0]  # Use first K-point as reference
        
        # Shift moiré BZ to k_ref
        moire_shifted = moire_vertices + k_ref
        moire_shifted_closed = np.vstack([moire_shifted, moire_shifted[0]])
        ax3.plot(moire_shifted_closed[:, 0], moire_shifted_closed[:, 1], 
                'green', linewidth=2, label='Moiré BZ @ k_ref')
        ax3.fill(moire_shifted_closed[:, 0], moire_shifted_closed[:, 1], 
                alpha=0.2, color='lightgreen')
        
        # Plot k_ref point
        ax3.scatter(k_ref[0], k_ref[1], c='red', s=120, marker='*', 
                   zorder=6, label='k_ref', edgecolors='black', linewidth=1)
        
        # Plot rcell vectors from k_ref
        if moire_rcell is not None:
            ax3.arrow(k_ref[0], k_ref[1], moire_rcell[0, 0], moire_rcell[0, 1], 
                     head_width=0.002, head_length=0.002, fc='red', ec='red', alpha=0.7)
            ax3.arrow(k_ref[0], k_ref[1], moire_rcell[1, 0], moire_rcell[1, 1], 
                     head_width=0.002, head_length=0.002, fc='orange', ec='orange', alpha=0.7)
        
        # Also show G-vectors around k_ref
        try:
            gvecs = np.loadtxt('g_vectors_debug.dat')
            ax3.scatter(gvecs[:, 0], gvecs[:, 1], c='purple', s=15, alpha=0.7, 
                       label=f'G-vectors', zorder=4)
        except:
            pass
    
    ax3.set_aspect('equal')
    ax3.grid(True, alpha=0.3)
    ax3.legend(fontsize=8)
    ax3.set_xlabel('kx (1/Å)')
    ax3.set_ylabel('ky (1/Å)')
    
    # ========== SUBPLOT 4: K-path in Graphene BZ ==========
    ax4 = fig.add_subplot(gs[1, 0])
    ax4.set_title('K-path in Graphene BZ', fontweight='bold')
    
    # Plot graphene BZ again (first layer only)
    bz_vertices, k_point, b1, b2 = parse_bz_file(filenames[0])
    if bz_vertices is not None:
        bz_closed = np.vstack([bz_vertices, bz_vertices[0]])
        ax4.plot(bz_closed[:, 0], bz_closed[:, 1], 'blue', linewidth=2, label='Graphene BZ')
        ax4.fill(bz_closed[:, 0], bz_closed[:, 1], alpha=0.15, color='lightblue')
    
    # Load and plot k-path
    kpath_frac, kpath_abs = parse_kpath_data('kpath_debug')
    if kpath_abs is not None:
        ax4.plot(kpath_abs[:, 0], kpath_abs[:, 1], 'ro-', linewidth=2, markersize=6, 
                label='K-path (absolute)')
        
        # Label special points
        if len(kpath_abs) >= 4:  # Based on your 4-point path
            labels = ['K', 'Γ', 'M?', "K'"]
            for i, (point, label) in enumerate(zip(kpath_abs, labels)):
                ax4.annotate(label, (point[0], point[1]), xytext=(5, 5), 
                           textcoords='offset points', fontsize=10, fontweight='bold')
    
    ax4.scatter(0, 0, c='black', s=80, marker='o', zorder=5, label='Origin')
    ax4.set_aspect('equal')
    ax4.grid(True, alpha=0.3)
    ax4.legend(fontsize=8)
    ax4.set_xlabel('kx (1/Å)')
    ax4.set_ylabel('ky (1/Å)')
    
    # ========== SUBPLOT 5: K-path in Moiré BZ ==========
    ax5 = fig.add_subplot(gs[1, 1])
    ax5.set_title('K-path in Moiré BZ', fontweight='bold')
    
    # Plot moiré BZ at origin
    if moire_vertices is not None:
        moire_closed = np.vstack([moire_vertices, moire_vertices[0]])
        ax5.plot(moire_closed[:, 0], moire_closed[:, 1], 'green', linewidth=2, label='Moiré BZ')
        ax5.fill(moire_closed[:, 0], moire_closed[:, 1], alpha=0.2, color='lightgreen')
    
    # Plot k-path (should be much smaller scale)
    if kpath_abs is not None:
        # The k-path in moiré coordinates (fractional * rcell)
        ax5.plot(kpath_abs[:, 0], kpath_abs[:, 1], 'ro-', linewidth=2, markersize=6, 
                label='K-path', alpha=0.8)
        
        # This might be very small compared to graphene BZ, so let's also show zoomed version
        ax5.set_xlim(-0.1, 0.1)  # Zoom in to see moiré-scale features
        ax5.set_ylim(-0.1, 0.1)
    
    ax5.scatter(0, 0, c='black', s=80, marker='o', zorder=5, label='Origin')
    ax5.set_aspect('equal')
    ax5.grid(True, alpha=0.3)
    ax5.legend(fontsize=8)
    ax5.set_xlabel('kx (1/Å)')
    ax5.set_ylabel('ky (1/Å)')
    
    # ========== SUBPLOT 6: Fractional vs Absolute Coordinates ==========
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.set_title('Fractional vs Absolute K-path', fontweight='bold')
    
    if kpath_frac is not None and kpath_abs is not None:
        # Plot fractional coordinates
        ax6_frac = ax6
        ax6_frac.plot(kpath_frac[:, 0], kpath_frac[:, 1], 'bo-', linewidth=2, 
                     markersize=6, label='Fractional coords')
        
        # Create second y-axis for absolute coordinates (different scale)
        ax6_abs = ax6.twinx()
        ax6_abs.plot(range(len(kpath_abs)), np.linalg.norm(kpath_abs, axis=1), 
                    'ro-', linewidth=2, markersize=6, label='|k| absolute')
        
        ax6_frac.set_xlabel('Point index')
        ax6_frac.set_ylabel('Fractional coordinates', color='blue')
        ax6_abs.set_ylabel('|k| (1/Å)', color='red')
        ax6_frac.grid(True, alpha=0.3)
        
        # Add point labels
        labels = ['K', 'Γ', 'M?', "K'"]
        for i, label in enumerate(labels[:len(kpath_frac)]):
            ax6_frac.annotate(label, (i, kpath_frac[i, 1]), xytext=(5, 5), 
                            textcoords='offset points', fontsize=10)
    
    plt.suptitle('TAPW Debugging: 2-Layer System', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.show()
    
    # Print diagnostic information
    print("\n" + "="*60)
    print("TAPW DEBUGGING SUMMARY")
    print("="*60)
    
    if len(all_k_points) > 0:
        print(f"K-points found: {len(all_k_points)}")
        for i, k_pt in enumerate(all_k_points):
            print(f"  Layer {i+1} K-point: ({k_pt[0]:.6f}, {k_pt[1]:.6f}), |K| = {np.linalg.norm(k_pt):.6f}")
    
    if kpath_abs is not None:
        print(f"\nK-path points: {len(kpath_abs)}")
        labels = ['K', 'Γ', 'M?', "K'"]
        for i, (point, label) in enumerate(zip(kpath_abs, labels)):
            print(f"  {label}: ({point[0]:.6f}, {point[1]:.6f}), |k| = {np.linalg.norm(point):.6f}")
    
    if moire_rcell is not None:
        print(f"\nMoiré rcell vectors:")
        print(f"  rcell[1]: ({moire_rcell[0, 0]:.6f}, {moire_rcell[0, 1]:.6f}), |r1| = {np.linalg.norm(moire_rcell[0]):.6f}")
        print(f"  rcell[2]: ({moire_rcell[1, 0]:.6f}, {moire_rcell[1, 1]:.6f}), |r2| = {np.linalg.norm(moire_rcell[1]):.6f}")
        
        # Check angle between moiré vectors
        dot_product = np.dot(moire_rcell[0], moire_rcell[1])
        magnitudes = np.linalg.norm(moire_rcell[0]) * np.linalg.norm(moire_rcell[1])
        if magnitudes > 0:
            angle = np.arccos(np.clip(dot_product / magnitudes, -1, 1)) * 180 / np.pi
            print(f"  Angle between rcell vectors: {angle:.1f}°")
    
    try:
        gvecs = np.loadtxt('g_vectors_debug.dat')
        print(f"\nG-vectors: {len(gvecs)} total")
        print(f"  Range: x=[{gvecs[:, 0].min():.3f}, {gvecs[:, 0].max():.3f}], y=[{gvecs[:, 1].min():.3f}, {gvecs[:, 1].max():.3f}]")
        print(f"  First 3 G-vectors:")
        for i in range(min(3, len(gvecs))):
            print(f"    G({i+1}): ({gvecs[i, 0]:.6f}, {gvecs[i, 1]:.6f})")
    except:
        print("\nNo G-vectors loaded")
    
    return fig

# Usage
if __name__ == "__main__":
    print("Starting comprehensive TAPW debugging visualization...")
    fig = plot_comprehensive_tapw_debug()
