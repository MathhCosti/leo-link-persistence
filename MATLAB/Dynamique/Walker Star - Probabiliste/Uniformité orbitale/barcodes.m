clear; clc; close all;

%% ============================================================
% BARCODE ZIGZAG H0 - WALKER STAR EN DEUX PARTIES
%
% Entree :
%   leo_zigzag_analysis_results.mat
%
% Sortie :
%   leo_H0_zigzag_barcodes_walker_star_deux_parts.mat
%% ============================================================

analysis_file = 'leo_zigzag_analysis_results.mat';
assert(isfile(analysis_file), 'Fichier introuvable : %s', analysis_file);

S = load(analysis_file, 'ZigzagAdjacency', 'ZigzagLabels', 'time_values');
ZigzagAdjacency = S.ZigzagAdjacency;
ZigzagLabels = S.ZigzagLabels;
time_values = S.time_values(:);

Nz = numel(ZigzagAdjacency);
fprintf('Nombre d''objets dans le zigzag : %d\n', Nz);

%% Conversion des indices zigzag en temps physiques
ZigzagTime = zeros(Nz,1);
for k = 1:Nz
    lab = ZigzagLabels(k);
    if abs(lab-round(lab)) < 1e-12
        idx = round(lab);
        ZigzagTime(k) = time_values(idx);
    else
        idx = floor(lab);
        ZigzagTime(k) = 0.5*(time_values(idx)+time_values(idx+1));
    end
end

%% Espaces H0
component_labels = cell(Nz,1);
h0_dims = zeros(Nz,1);
for k = 1:Nz
    comp = conncomp(graph(ZigzagAdjacency{k}));
    comp = comp(:);
    component_labels{k} = comp;
    h0_dims(k) = max(comp);
end
fprintf('Dimensions H0 min/max : %d / %d\n', min(h0_dims), max(h0_dims));

%% Applications du module zigzag
maps = cell(Nz-1,1);
for k = 1:Nz-1
    if mod(k,2)==1
        maps{k}.type = 'f';
        maps{k}.mat = build_H0_map(component_labels{k}, component_labels{k+1}, ...
            h0_dims(k), h0_dims(k+1));
    else
        maps{k}.type = 'g';
        maps{k}.mat = build_H0_map(component_labels{k+1}, component_labels{k}, ...
            h0_dims(k+1), h0_dims(k));
    end
end

%% Barcode zigzag H0
intervals_all = zigzag_barcode_from_module_mod2(h0_dims, maps);
birth_index_all = intervals_all(:,1);
death_index_all = intervals_all(:,2);
birth_time_all = ZigzagTime(birth_index_all);
death_time_all = ZigzagTime(death_index_all);
lifetimes_all = death_time_all-birth_time_all;

% Retrait de la composante globale persistante pour les statistiques et p_disp.
persistent_idx = (birth_index_all==1) & (death_index_all==Nz);
fprintf('Barres globales persistantes retirees : %d\n', nnz(persistent_idx));
keep = ~persistent_idx;

intervals = intervals_all(keep,:);
birth_index = birth_index_all(keep);
death_index = death_index_all(keep);
birth_time = birth_time_all(keep);
death_time = death_time_all(keep);
lifetimes = lifetimes_all(keep);

output_file = 'leo_H0_zigzag_barcodes_walker_star_deux_parts.mat';
save(output_file, 'intervals','birth_index','death_index', ...
    'birth_time','death_time','lifetimes','ZigzagTime','ZigzagLabels', ...
    'h0_dims','intervals_all','birth_index_all','death_index_all', ...
    'birth_time_all','death_time_all','lifetimes_all','persistent_idx');
fprintf('Barcodes sauvegardes dans %s\n', output_file);

%% Affichage des plus longues barres
[~,order] = sort(lifetimes,'descend');
maxBarsToPlot = 150;
order = order(1:min(maxBarsToPlot,numel(order)));
figure; hold on; grid on;
for ii = 1:numel(order)
    id = order(ii);
    if abs(death_time(id)-birth_time(id)) < 1e-12
        plot(birth_time(id),ii,'ko','MarkerSize',4);
    else
        plot([birth_time(id),death_time(id)],[ii,ii],'k-','LineWidth',1.2);
    end
end
xlabel('Temps (s)'); ylabel('Barres H_0 triees par duree decroissante');
title(sprintf('Barcode zigzag H_0 - Walker Star deux parties - %d barres',numel(order)));
hold off;

figure;
histogram(lifetimes,30); grid on;
xlabel('Duree de vie (s)'); ylabel('Nombre de barres');
title('Distribution des durees de vie H_0 - Walker Star deux parties');

%% ============================================================
%  PROBABILITE EMPIRIQUE DE SURVIE DES BARRES H0
%
%  S(L) = P(T >= L)
%
%  Les barres qui atteignent la fin de la fenetre d'observation sont
%  considerees comme censurees a droite et ne sont pas utilisees pour
%  estimer directement la distribution des durees terminees.
%% ============================================================

T_end = max(ZigzagTime);
tol_survival = 1e-10 * max(1, abs(T_end));

is_right_censored = abs(death_time - T_end) <= tol_survival;
completed_lifetimes = lifetimes(~is_right_censored & lifetimes > 0);

if isempty(completed_lifetimes)
    warning(['Aucune duree positive completement observee : ', ...
             'la courbe de survie ne peut pas etre calculee.']);
    survival_times = [];
    survival_emp = [];
else
    survival_times = unique(sort(completed_lifetimes));

    survival_emp = zeros(size(survival_times));
    for i = 1:numel(survival_times)
        survival_emp(i) = mean(completed_lifetimes >= survival_times(i));
    end

    figure;
    stairs([0; survival_times], [1; survival_emp], ...
        'LineWidth', 1.8);
    grid on;
    xlabel('Duree L (s)');
    ylabel('Probabilite de survie S(L) = P(T \geq L)');
    title('Probabilite empirique de duree de vie des barres H_0');
    ylim([0 1.05]);

    figure;
    semilogy(survival_times, survival_emp, 'o-', ...
        'LineWidth', 1.5, 'MarkerSize', 4);
    grid on;
    xlabel('Duree L (s)');
    ylabel('Probabilite de survie S(L)');
    title('Survie empirique des barres H_0 - echelle semi-logarithmique');
end

fprintf('\nAnalyse de survie H0 :\n');
fprintf('Barres censurees a droite : %d\n', nnz(is_right_censored));
fprintf('Durees positives terminees : %d\n', numel(completed_lifetimes));

save(output_file, ...
    'T_end', 'is_right_censored', 'completed_lifetimes', ...
    'survival_times', 'survival_emp', ...
    '-append');

fprintf('\nStatistiques H0 hors composante globale :\n');
if isempty(lifetimes)
    fprintf('Aucune barre non persistante.\n');
else
    fprintf('Nombre de barres : %d\n',numel(lifetimes));
    fprintf('Duree moyenne    : %.2f s\n',mean(lifetimes));
    fprintf('Duree mediane    : %.2f s\n',median(lifetimes));
    fprintf('Duree maximale   : %.2f s\n',max(lifetimes));
end

%% Fonctions locales
function M = build_H0_map(labels_source,labels_target,dim_source,dim_target)
M = zeros(dim_target,dim_source);
for c = 1:dim_source
    vertices = find(labels_source==c);
    target_comps = unique(labels_target(vertices));
    if numel(target_comps)~=1
        error('Inclusion invalide en H0 : une composante source a plusieurs images.');
    end
    M(target_comps(1),c)=1;
end
M = mod(M,2);
end

function intervals = zigzag_barcode_from_module_mod2(dims,maps)
n = numel(dims);
R = {zeros(dims(1),0); eye(dims(1))};
b = 1;
r = filtration_quotient_dims(R);
intervals = zeros(0,2);
for k = 1:n-1
    if maps{k}.type=='f'
        M = maps{k}.mat;
        Rnext = cell(numel(R)+1,1);
        for i=1:numel(R), Rnext{i}=gf2_col_basis(M*R{i}); end
        Rnext{end}=eye(dims(k+1));
        bnext=[b,k+1];
        rnext=filtration_quotient_dims(Rnext);
        for i=1:numel(r)
            c=r(i)-rnext(i);
            if c>0, intervals=[intervals;repmat([b(i),k],c,1)]; end %#ok<AGROW>
        end
    else
        N=maps{k}.mat;
        Rnext=cell(numel(R)+1,1);
        Rnext{1}=zeros(dims(k+1),0);
        for i=1:numel(R), Rnext{i+1}=gf2_preimage(N,R{i}); end
        bnext=[k+1,b];
        rnext=filtration_quotient_dims(Rnext);
        for i=1:numel(r)
            c=r(i)-rnext(i+1);
            if c>0, intervals=[intervals;repmat([b(i),k],c,1)]; end %#ok<AGROW>
        end
    end
    R=Rnext; b=bnext; r=rnext;
end
for i=1:numel(r)
    if r(i)>0, intervals=[intervals;repmat([b(i),n],r(i),1)]; end %#ok<AGROW>
end
end

function dims = filtration_quotient_dims(R)
dims=zeros(1,numel(R)-1);
for i=1:numel(dims), dims(i)=gf2_rank(R{i+1})-gf2_rank(R{i}); end
end

function P = gf2_preimage(A,S)
A=mod(full(A),2); S=mod(full(S),2); n=size(A,2);
Z=gf2_null([A,S]); P=gf2_col_basis(Z(1:n,:));
end

function B = gf2_col_basis(A)
A=mod(full(A),2);
if isempty(A), B=zeros(size(A,1),0); return; end
[~,p]=gf2_rref(A);
if isempty(p), B=zeros(size(A,1),0); else, B=A(:,p); end
end

function r = gf2_rank(A)
if isempty(A), r=0; return; end
[~,p]=gf2_rref(mod(full(A),2)); r=numel(p);
end

function Z = gf2_null(A)
A=mod(full(A),2); [R,pivots]=gf2_rref(A); n=size(A,2);
free_cols=setdiff(1:n,pivots);
Z=zeros(n,numel(free_cols));
for j=1:numel(free_cols)
    f=free_cols(j); z=zeros(n,1); z(f)=1;
    for p=1:numel(pivots), z(pivots(p))=R(p,f); end
    Z(:,j)=mod(z,2);
end
end

function [R,pivots] = gf2_rref(A)
R=mod(full(A),2); [m,n]=size(R); pivots=[]; row=1;
for col=1:n
    if row>m, break; end
    q=find(R(row:m,col),1);
    if isempty(q), continue; end
    q=q+row-1;
    if q~=row, tmp=R(row,:); R(row,:)=R(q,:); R(q,:)=tmp; end
    for rr=1:m
        if rr~=row && R(rr,col)==1, R(rr,:)=mod(R(rr,:)+R(row,:),2); end
    end
    pivots(end+1)=col; %#ok<AGROW>
    row=row+1;
end
end
