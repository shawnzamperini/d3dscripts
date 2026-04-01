# Provided by Tess on 4/1/2026
# Runs from command line

import numpy as np
import postgkyl as pg
import fnmatch
import os
import re

# --- 1. Exact p=1 3D Serendipity Basis Evaluator ---
def eval_p1_basis(xi, eta, zeta):
    """
    Evaluates the 8 orthonormal basis functions for p=1 serendipity in 3D 
    at logical coordinates (xi, eta, zeta) in [-1, 1].
    """
    return np.array([
        1.0 / np.sqrt(8.0),
        np.sqrt(3.0/8.0) * xi,
        np.sqrt(3.0/8.0) * eta,
        np.sqrt(3.0/8.0) * zeta,
        np.sqrt(9.0/8.0) * xi * eta,
        np.sqrt(9.0/8.0) * xi * zeta,
        np.sqrt(9.0/8.0) * eta * zeta,
        np.sqrt(27.0/8.0) * xi * eta * zeta
    ])

# --- 2. Extract Raw Modal Data at Specific Nodes ---
def get_raw_nodal_values(filename, comp, local_nodes):
    """
    Loads raw DG modal coefficients and evaluates them exactly at the 
    specified local quadrature nodes for a specific vector component.
    """
    gdata = pg.data.GData(filename)
    modal_data = gdata.get_values() # Shape: (Nx, Ny, Nz, 72)
    
    # Extract the 8 basis coefficients for the specific physical component 'comp'
    coeffs = modal_data[..., comp*8 : (comp+1)*8]
    
    # Evaluate at all provided local nodes
    nodal_values = []
    for (xi, eta, zeta) in local_nodes:
        basis_vals = eval_p1_basis(xi, eta, zeta)
        # Dot product of modal coefficients with basis functions
        val_at_node = np.sum(coeffs * basis_vals, axis=-1)
        nodal_values.append(val_at_node)
        
    return np.stack(nodal_values, axis=-1) # Shape: (Nx, Ny, Nz, num_nodes)

def find_prefix(pattern, path):
    for name in os.listdir(path):
        if fnmatch.fnmatch(name, '*' + pattern):
            return re.sub(pattern, '', name)
    raise FileNotFoundError("ERROR: file prefix not found!")

# --- 3. Main Matrix Orthogonality Check ---
def check_orthogonality_matrix_nodes(file_prefix):
    # Quadrature nodes (e.g., Gauss-Lobatto corners for standard p=1 mappings)
    c = 1.0/np.sqrt(3.0) # For Gauss-Legendre; use 1.0 for Gauss-Lobatto
    local_quad_nodes = [
        (-c, -c, -c), ( c, -c, -c), (-c,  c, -c), ( c,  c, -c),
        (-c, -c,  c), ( c, -c,  c), (-c,  c,  c), ( c,  c,  c)
    ]
    
    print(f"Loading 9 components and evaluating at {len(local_quad_nodes)} nodes...")
    
    dxdz_comps = []
    dzdx_comps = []
    
    # Loop to load all 9 components
    for i in range(9):
        comp_dxdz = get_raw_nodal_values(f"{file_prefix}-dxdz.gkyl", i, local_quad_nodes)
        comp_dzdx = get_raw_nodal_values(f"{file_prefix}-dzdx.gkyl", i, local_quad_nodes)
        
        dxdz_comps.append(comp_dxdz)
        dzdx_comps.append(comp_dzdx)

    print("Stacking and reshaping into 3x3 matrices...")
    # Stack along a new last axis. Shape: (Nx, Ny, Nz, num_nodes, 9)
    dxdz_flat = np.stack(dxdz_comps, axis=-1)
    dzdx_flat = np.stack(dzdx_comps, axis=-1)
    
    # Reshape the 9 components into a 3x3 matrix at each spatial node
    # Shape: (Nx, Ny, Nz, num_nodes, 3, 3)
    J = dxdz_flat.reshape(*dxdz_flat.shape[:-1], 3, 3)
    J_inv = dzdx_flat.reshape(*dzdx_flat.shape[:-1], 3, 3)
    
    print("Computing J @ J_inv^T across all nodes...")
    # Transpose the inner 3x3 matrices of J_inv to dot rows with rows
    identity_check = J @ J_inv.swapaxes(-1, -2)
    
    # Create 3x3 identity matrix for comparison
    I_ideal = np.eye(3)
    
    # Calculate the absolute error everywhere
    error_matrix = np.abs(identity_check - I_ideal)
    max_error = np.max(error_matrix)
    
    # Find the specific index and node of the maximum error
    max_err_idx = np.unravel_index(np.argmax(error_matrix), error_matrix.shape)
    ix, iy, iz, node_idx, row, col = max_err_idx
    
    print("-" * 60)
    print(f"Maximum orthogonality error at exact nodes: {max_error:.2e}")
    print(f"Grid Index:           (Nx={ix}, Ny={iy}, Nz={iz})")
    print(f"Local Node Index:     {node_idx} {local_quad_nodes[node_idx]}")
    print(f"Specific Matrix Error: J[{row}, :] dot J_inv[{col}, :]")
    print("-" * 60)
    
    if max_error < 1e-10:
        print("Success! Your dual and tangent vectors are perfectly orthogonal everywhere.")
    else:
        print("Warning: Error is larger than machine precision. If isolated, check for grid singularities.")
        
    return max_error

# --- 4. Execute ---
if __name__ == "__main__":
    my_file_prefix = find_prefix('-dxdz.gkyl', '.')
    check_orthogonality_matrix_nodes(my_file_prefix)
