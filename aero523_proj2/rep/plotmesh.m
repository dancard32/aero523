function plotmesh

Mesh = readgri('mesh0.gri');

figure(1); clf;
for e = 1:Mesh.nElem,
  I = Mesh.Elem(e,:); I = [I,I(1)];
  plot(Mesh.Node(I,1), Mesh.Node(I,2), 'k-'); hold on;
end
axis equal; axis off
