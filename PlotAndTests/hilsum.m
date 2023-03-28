n=logspace(16,21,101);
nv=[1e17 5e18 ];

for k=1:length(nv)
n0=nv(k);
y(k,:)=1./(1+(n./n0).^.35);
yn(k,:)=1./(1+(n./n0).^.5);
end

figure,
loglog(n,y,'linewidth',2), hold on,
ax = gca;
ax.ColorOrderIndex = 1;
loglog(n,yn,'--','linewidth',2),
ylim([.1 1])
grid

pausak

figure, semilogx(n,y), hold on,
ax = gca;
ax.ColorOrderIndex = 1;
semilogx(n,yn,'--'),
