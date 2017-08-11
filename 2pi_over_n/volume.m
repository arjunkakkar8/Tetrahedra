% Assume that dihedral angles and face areas are already computed
% in another Matlab file

surfaceareas = abs(sum(areas, 2));

surfaceareas

cayleymenger = @(x) [
    0,      x(1)^2, x(2)^2, x(3)^2, 1;
    x(1)^2, 0,      x(4)^2, x(5)^2, 1;
    x(2)^2, x(4)^2, 0,      x(6)^2, 1;
    x(3)^2, x(5)^2, x(6)^2, 0,      1;
    1,      1,      1,      1,      0];

vol = @(x) sqrt(det(cayleymenger(x)) / 288);

% Computed by Python
edges = vpa([
0.7158639394115428, 1.3374806099528422, 0.7861513777574243, 1.2720196495140668, 1.1582921852882702, 1.8741561246998029;
1.1582921852882686, 0.7861513777574232, 0.7861513777574234, 1.2720196495140685, 1.272019649514069, 1.3374806099528451;
1.0040942904713324, 1.1026817617063975, 1.1026817617063966, 0.681494807509166, 1.6246586898923236, 1.875994263901196;
0.974003746425297, 0.9740037464252967, 0.9740037464252967, 1.1246826503806981, 1.1246826503806981, 1.5905414575341017;
1.3433178570325923, 0.8302160933408144, 0.9117312518252875, 0.9117312518252881, 1.4752121540588046, 1.3433178570325937;
0.5197025656447412, 0.8408964152537145, 1.360598980898456, 0.9885330231834818, 1.456475315121971, 1.5994800305125612;
0.7071067811865475, 0.9999999999999998, 0.9999999999999998, 1.224744871391589, 1.224744871391589, 1.4142135623730954;
0.9036020036098451, 1.2778862084925453, 0.9036020036098447, 0.903602003609845, 1.2778862084925446, 1.5650845800732873;
0.9114860908270551, 0.866874786193467, 1.2772276674444207, 0.866874786193467, 1.2772276674444207, 1.4026328680513274;
1.2194228708043555, 1.1597400673978266, 0.6526755256162672, 1.1597400673978249, 1.159740067397825, 1.2194228708043562;
1.189207115002721, 1.0298835719535588, 1.0298835719535588, 1.0298835719535588, 1.0298835719535588, 1.189207115002721
]);

volumes = [];

% Compute volumes
for i=1:length(edges)
    volumes = [volumes; vol(edges(i,:))];
end

volumes

% Compute normalized surface areas
normalizedsurfaceareas = surfaceareas ./ (volumes .^ (2/3))