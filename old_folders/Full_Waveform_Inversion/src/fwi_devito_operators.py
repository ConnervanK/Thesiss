import numpy as np
from devito import TimeFunction, Function, Eq, Operator, solve, configuration
import devito

configuration['log-level'] = 'INFO'

def setup_devito_fwi(shape=(500, 250), origin=(0., 0.), spacing=(0.002, 0.002),
                     vp=0.1224, space_order=4, time_order=2):
    """
    Sets up the Devito FWI operators.
    shape: (nx, ny) grid
    spacing: (dx, dy) in meters
    vp: initial velocity in m/ns
    """
    # Create the Devito Grid
    from devito import Grid
    grid = Grid(shape=shape, extent=(shape[0]*spacing[0], shape[1]*spacing[1]), origin=origin)
    
    # Velocity field (m/ns)
    m = Function(name='m', grid=grid)
    # Background model (1 / vp^2)
    m.data[:] = 1.0 / (vp**2)
    
    # Forward wavefield
    u = TimeFunction(name='u', grid=grid, time_order=time_order, space_order=space_order)
    
    # Adjoint wavefield
    v = TimeFunction(name='v', grid=grid, time_order=time_order, space_order=space_order)
    
    # The acoustic wave equation in Devito is: m * d^2 u / dt^2 - laplacian(u) = src
    # And adjoint: m * d^2 v / dt^2 - laplacian(v) = data_residual
    return grid, m, u, v

def forward_operator(grid, m, u, src, rec, space_order=4):
    """
    Returns the forward operator.
    src, rec are Devito PointSource/Receiver instances.
    """
    # Acoustic wave equation
    pde = m * u.dt2 - u.laplace
    # Dampening/absorber boundaries usually applied, but let's assume simplest BCs to start,
    # or rely on Devito's default damping if created with boundary layers. 
    # Let's keep it simple.
    stencil = Eq(u.forward, solve(pde, u.forward))
    
    # Source injection
    src_term = src.inject(field=u.forward, expr=src * grid.time_dim.spacing**2 / m)
    
    # Receiver interpolation
    rec_term = rec.interpolate(expr=u.forward)
    
    op = Operator([stencil] + src_term + rec_term, subs=grid.spacing_map, name='Forward')
    return op

def adjoint_operator(grid, m, v, srca, rec, space_order=4):
    """
    Adjoint equation solves backwards in time from residual injected at receivers.
    """
    pde = m * v.dt2 - v.laplace
    # Backward in time!
    stencil = Eq(v.backward, solve(pde, v.backward))
    
    # Inject data residual at receivers
    rec_term = rec.inject(field=v.backward, expr=rec * grid.time_dim.spacing**2 / m)
    
    # Extract at source (optional, if we were doing source estimation)
    srca_term = srca.interpolate(expr=v.backward) if srca else []
    
    op = Operator([stencil] + rec_term + srca_term, subs=grid.spacing_map, name='Adjoint')
    return op

def gradient_operator(grid, u, v, m):
    """
    The gradient of the objective function w.r.t model m.
    grad = sum_t -(v.dt2 * u)
    """
    grad = Function(name='grad', grid=grid)
    # The gradient operator is usually executed within loops or as a separate update
    gradient_eq = Eq(grad, grad - u.dt2 * v)
    op = Operator([gradient_eq], name='Gradient')
    return op, grad
