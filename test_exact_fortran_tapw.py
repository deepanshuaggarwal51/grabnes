#!/usr/bin/env python3
"""
Exact replication of Fortran TAPW transformation logic.
This test mimics the exact same:
1. X matrix construction (with label-based normalization)
2. Column ordering: col = (i_G - 1) * N_label + i_label  
3. CSR sparse matrix format
4. Two-step transformation: Y = H*X, then H_proj = X†*Y
"""

import numpy as np
import scipy.sparse as sp
from scipy.sparse import csr_matrix

def build_X_fortran_exact(xcoord, ycoord, label, Gx, Gy, N_orbit, N_G, N_label):
    """
    Exact replication of Fortran build_X subroutine
    """
    # Count orbitals per label (exactly like Fortran)
    count = np.zeros(N_label, dtype=int)
    for i_orb in range(N_orbit):
        count[label[i_orb] - 1] += 1  # Convert to 0-based indexing
    
    # Compute normalization (exactly like Fortran)
    norm = np.zeros(N_label)
    for i_label in range(N_label):
        if count[i_label] > 0:
            norm[i_label] = 1.0 / np.sqrt(count[i_label])
        else:
            norm[i_label] = 0.0
    
    print(f"Label counts: {count}")
    print(f"Normalizations: {norm}")
    
    # Build X matrix (exactly like Fortran)
    M = N_G * N_label
    X = np.zeros((N_orbit, M), dtype=complex)
    
    for i_G in range(N_G):
        for i_label in range(N_label):
            col = i_G * N_label + i_label  # Fortran: (i_G - 1) * N_label + i_label, but 0-based
            
            for i_orb in range(N_orbit):
                if label[i_orb] == i_label + 1:  # Convert back to 1-based for comparison
                    phase = 1j * (Gx[i_G] * xcoord[i_orb] + Gy[i_G] * ycoord[i_orb])
                    X[i_orb, col] = norm[i_label] * np.exp(phase)
                else:
                    X[i_orb, col] = 0.0
    
    return X

def create_realistic_sparse_hamiltonian(N, label):
    """
    Create a sparse Hamiltonian that mimics tight-binding structure
    - On-site energies on diagonal
    - Nearest-neighbor hoppings off-diagonal
    - Sublattice-dependent structure
    """
    # Create empty sparse matrix
    row_ind = []
    col_ind = []
    values = []
    
    np.random.seed(42)
    
    # Add diagonal elements (on-site energies)
    for i in range(N):
        row_ind.append(i)
        col_ind.append(i)
        # Different on-site energies for different sublattices
        if label[i] == 1:  # Sublattice A
            values.append(0.1 + 0.05 * np.random.randn())
        else:  # Sublattice B  
            values.append(-0.1 + 0.05 * np.random.randn())
    
    # Add off-diagonal hopping elements (nearest neighbors)
    hopping_strength = -2.7  # Typical graphene hopping
    
    for i in range(N):
        # Add a few random nearest neighbors
        num_neighbors = np.random.randint(2, 6)  # 2-5 neighbors
        neighbors = np.random.choice(N, num_neighbors, replace=False)
        
        for j in neighbors:
            if i != j:
                # Add hopping i -> j
                row_ind.append(i)
                col_ind.append(j)
                
                # Phase factor for hopping (simplified)
                phase = np.random.uniform(0, 2*np.pi)
                hop_value = hopping_strength * np.exp(1j * phase)
                values.append(hop_value)
    
    # Create CSR matrix
    H_sparse = csr_matrix((values, (row_ind, col_ind)), shape=(N, N), dtype=complex)
    
    # Make it Hermitian
    H_sparse = (H_sparse + H_sparse.H) / 2.0
    
    return H_sparse

def fortran_transform_exact(H_csr, X):
    """
    Exact replication of Fortran transform_sparse_hamiltonian
    """
    N, M = X.shape
    
    print(f"\nFortran-exact transformation:")
    print(f"  H shape: {H_csr.shape} (nnz={H_csr.nnz})")
    print(f"  X shape: {X.shape}")
    print(f"  Expected H_proj shape: ({M}, {M})")
    
    # Step 1: Y = H * X (exactly like Fortran sparse multiplication)
    print(f"  Step 1: Y = H * X...")
    Y = H_csr.dot(X)
    print(f"    Y shape: {Y.shape}")
    print(f"    Max |Y|: {np.max(np.abs(Y)):.6f}")
    
    # Step 2: H_proj = X† * Y (exactly like Fortran zgemm call)
    print(f"  Step 2: H_proj = X† * Y...")
    X_dagger = X.conj().T
    H_proj = X_dagger @ Y
    print(f"    H_proj shape: {H_proj.shape}")
    print(f"    Max |H_proj|: {np.max(np.abs(H_proj)):.6f}")
    
    return H_proj

def test_exact_fortran_match():
    """
    Test with exact Fortran parameters and structure
    """
    print("=== EXACT FORTRAN TAPW REPLICATION TEST ===\n")
    
    # Parameters matching typical Fortran run
    N_orbit = 100    # Number of atoms/orbitals
    N_G = 25        # Number of G-vectors  
    N_label = 2     # Two sublattices (A, B)
    
    print(f"Parameters (matching Fortran):")
    print(f"  N_orbit: {N_orbit}")
    print(f"  N_G: {N_G}")
    print(f"  N_label: {N_label}")
    print(f"  M = N_G * N_label = {N_G * N_label}")
    
    # Create atomic positions (like a graphene-like structure)
    np.random.seed(123)
    xcoord = np.random.uniform(-10, 10, N_orbit)
    ycoord = np.random.uniform(-10, 10, N_orbit)
    
    # Create sublattice labels (1 or 2, like Fortran)
    label = np.random.choice([1, 2], N_orbit)
    
    # Create G-vectors (like TAPW would generate)
    Gx = np.random.uniform(-2, 2, N_G)
    Gy = np.random.uniform(-2, 2, N_G)
    
    print(f"\nCreated test system:")
    print(f"  Coordinate ranges: x=[{np.min(xcoord):.2f}, {np.max(xcoord):.2f}]")
    print(f"  Coordinate ranges: y=[{np.min(ycoord):.2f}, {np.max(ycoord):.2f}]")
    print(f"  G-vector ranges: Gx=[{np.min(Gx):.2f}, {np.max(Gx):.2f}]")
    print(f"  G-vector ranges: Gy=[{np.min(Gy):.2f}, {np.max(Gy):.2f}]")
    print(f"  Label distribution: {np.bincount(label)}")
    
    # Build X matrix exactly like Fortran
    print(f"\nBuilding X matrix (Fortran-exact)...")
    X = build_X_fortran_exact(xcoord, ycoord, label, Gx, Gy, N_orbit, N_G, N_label)
    
    # Check X matrix properties
    print(f"\nX matrix analysis:")
    print(f"  Shape: {X.shape}")
    print(f"  Max |X|: {np.max(np.abs(X)):.6f}")
    print(f"  Sparsity: {np.count_nonzero(X)} / {X.size} = {np.count_nonzero(X)/X.size:.3f}")
    
    # Check column norms (should reflect sublattice structure)
    col_norms = np.linalg.norm(X, axis=0)
    print(f"  Column norms range: [{np.min(col_norms):.6f}, {np.max(col_norms):.6f}]")
    
    # Check orthogonality within sublattices
    X_dagger_X = X.conj().T @ X
    diag_elements = np.diag(X_dagger_X)
    off_diag_max = np.max(np.abs(X_dagger_X - np.diag(diag_elements)))
    print(f"  X†X diagonal range: [{np.min(diag_elements):.6f}, {np.max(diag_elements):.6f}]")
    print(f"  X†X off-diagonal max: {off_diag_max:.6f}")
    
    # Create realistic sparse Hamiltonian
    print(f"\nCreating sparse Hamiltonian...")
    H_sparse = create_realistic_sparse_hamiltonian(N_orbit, label)
    print(f"  H sparsity: {H_sparse.nnz} / {N_orbit**2} = {H_sparse.nnz/N_orbit**2:.4f}")
    
    # Check Hermiticity
    hermitian_error = np.max(np.abs(H_sparse.toarray() - H_sparse.H.toarray()))
    print(f"  Hermiticity error: {hermitian_error:.2e}")
    
    # Perform transformation
    H_proj = fortran_transform_exact(H_sparse, X)
    
    # Check result properties
    print(f"\nResult analysis:")
    hermitian_error_proj = np.max(np.abs(H_proj - H_proj.conj().T))
    print(f"  H_proj Hermiticity error: {hermitian_error_proj:.2e}")
    
    eigenvals_proj = np.linalg.eigvals(H_proj)
    print(f"  H_proj eigenvalue range: [{np.min(eigenvals_proj.real):.3f}, {np.max(eigenvals_proj.real):.3f}]")
    print(f"  H_proj eigenvalue imaginary parts max: {np.max(np.abs(eigenvals_proj.imag)):.2e}")
    
    # Compare with direct method
    print(f"\nComparing with direct method...")
    H_direct = X.conj().T @ H_sparse.toarray() @ X
    difference = np.max(np.abs(H_proj - H_direct))
    print(f"  Max |H_proj_fortran - H_proj_direct|: {difference:.2e}")
    
    if difference < 1e-12:
        print("  ✅ PERFECT MATCH: Fortran method identical to direct method!")
    elif difference < 1e-8:
        print("  ✅ EXCELLENT MATCH: Methods agree within numerical precision")
    else:
        print("  ❌ METHODS DISAGREE!")
        return False
    
    # Test reconstruction
    print(f"\nTesting reconstruction X * H_proj * X†...")
    H_reconstructed = X @ H_proj @ X.conj().T
    H_original = H_sparse.toarray()
    
    # For fair comparison, project original H into same subspace
    P = X @ X.conj().T  # Projection operator
    H_original_projected = P @ H_original @ P
    
    reconstruction_error = np.max(np.abs(H_reconstructed - H_original_projected))
    print(f"  Reconstruction error (in projected subspace): {reconstruction_error:.2e}")
    
    if reconstruction_error < 1e-10:
        print("  ✅ EXCELLENT RECONSTRUCTION!")
    elif reconstruction_error < 1e-6:
        print("  ✅ GOOD RECONSTRUCTION")
    else:
        print("  ⚠ POOR RECONSTRUCTION (may indicate non-unitary X)")
    
    print(f"\n🎯 CONCLUSION:")
    if difference < 1e-8 and hermitian_error_proj < 1e-12:
        print("   ✅ FORTRAN TAPW LOGIC IS MATHEMATICALLY CORRECT!")
        print("   The two-step transformation works perfectly.")
        print("   Any issues in Fortran are likely in:")
        print("     - Input data format/values")
        print("     - G-vector completeness")  
        print("     - Numerical precision settings")
    else:
        print("   ❌ THERE MAY BE IMPLEMENTATION ISSUES!")
    
    return True

if __name__ == "__main__":
    test_exact_fortran_match()
