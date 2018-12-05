function tight_ax(ax)
% Removes the annoying extra L-R margins around plots

t = ax.TightInset;
ax.LooseInset = [t(1) 0 t(3) 0] *1.05;

end