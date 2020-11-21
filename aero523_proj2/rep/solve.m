function [u, err, V, E, BE, IE, u0] = solve()
    mesh = readgri('mesh0.gri');
    V = mesh.Node; E = mesh.Elem; BE = mesh.B.nodes;
    [IE, BE] = edgehash(E);
    
    
end

