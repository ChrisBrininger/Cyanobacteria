%% Plot %%
% 22C
x = [0,6,24];
PelletUrea = [9.080, 9.364, 9.439];
SupUrea = [9.380, 9.391, 9.522];
ContUrea = [9.136, 9.149, 9.386];
Pellet = [8.815, 8.996, 9.316];
Sup = [9.371, 9.343, 9.408];
Cont = [8.856, 8.998, 9.331];

plot (x, PelletUrea, "r", x, SupUrea, "g", x, ContUrea, "b", x, Pellet, "c", x, Sup, "m", x, Cont, "k")
legend ('Pellet + Urea','Supernatant + Urea', 'Fresh Media + Urea','Pellet - Urea','Supernatant - Urea', 'Fresh Media - Urea','Location','southeast')
title ('Room Temperature')

%% Plot %%
% 30C
x = [0,6,24];
PelletUrea = [9.122, 9.337, 9.434];
SupUrea = [9.417, 9.421, 9.489];
ContUrea = [9.181, 9.199, 9.453];
Pellet = [8.821, 9.103, 9.405];
Sup = [9.375, 9.346, 9.45];
Cont = [8.877, 9.14, 9.425];

plot (x, PelletUrea, "r", x, SupUrea, "g", x, ContUrea, "b", x, Pellet, "c", x, Sup, "m", x, Cont, "k")
legend ('Pellet + Urea','Supernatant + Urea', 'Fresh Media + Urea','Pellet - Urea','Supernatant - Urea', 'Fresh Media - Urea','Location','southeast')
title ('30 C')
