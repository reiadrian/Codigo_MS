function [Fext] = f_Fflux(Fext,u,e_DatSet,e_VG)
% Determina las fuerzas debidado a las poropresiones en cada paso de tiempo
% y los suma al vector global de fuerzas externas

nSet = e_VG.nSet;
c_Fext = cell(nSet,1);
c_FilFza = cell(nSet,1);


for iSet = 1:nSet       
%     e_DatMatSet = e_DatSet(iSet).e_DatMat;
%     if e_DatMatSet.conshyp==14
%     Dtime = e_VG.Dtime;
%     e_DatElemSet = e_DatSet(iSet).e_DatElem;
%     m_DerCa_p = e_DatSet(iSet).m_DerCa_p;
%     m_DetJT_p = e_DatSet(iSet).m_DetJT_p;
%     
%     nElem = e_DatSet(iSet).nElem;
% 
%     dofpe = e_DatElemSet.dofpe;
%     dofpe_p = e_DatElemSet.dofpe_p;
%     wg = e_DatElemSet.wg;
%     nPG = e_DatElemSet.npg;
%     pos_p = e_DatElemSet.pos_p;
%     pos_p0 = e_DatElemSet.pos_p0;
% 
%     PermK = e_DatMatSet.m_PermK;

    % Grados de libertad y coordenadas de los nodos de los elementos del set
    e_DatElemSet = e_DatSet(iSet).e_DatElem;
    e_DatMatSet = e_DatSet(iSet).e_DatMat;
    nElem = e_DatSet(iSet).nElem;
    m_DofElem = e_DatSet(iSet).m_DofElem;
    dofElemSet = m_DofElem(:); 
    m_FilFza = dofElemSet';
    uElemSet  = reshape(u(dofElemSet),[],nElem);
    dofpe = e_DatElemSet.dofpe;
    m_Fext = zeros(dofpe,nElem);
    if e_DatMatSet.conshyp==14
        
        Dtime = e_VG.Dtime;
        
        m_DerCa_p = e_DatSet(iSet).m_DerCa_p;
        m_DetJT_p = e_DatSet(iSet).m_DetJT_p;

        dofpe_p = e_DatElemSet.dofpe_p;
        wg = e_DatElemSet.wg;
        nPG = e_DatElemSet.npg;
        pos_p = e_DatElemSet.pos_p;
        pos_p0 = e_DatElemSet.pos_p0;

        PermK = e_DatMatSet.m_PermK;
        for iElem = 1:nElem
            Hww = zeros(dofpe_p,dofpe_p);
            m_Dercae_p = m_DerCa_p(:,:,:,iElem);
            ue = uElemSet(:,iElem);
            ue_p = ue(pos_p);
            m_pesoPG_p = m_DetJT_p(:,iElem).*wg;
        for iPG = 1:nPG
            DerivN = m_Dercae_p(:,:,iPG);
            Hww = Hww + DerivN'*PermK*DerivN*m_pesoPG_p(iPG); 
        end %for(iPG)

        m_Fext(pos_p,iElem) =  Dtime*Hww*ue_p;
        m_Fext(pos_p0,iElem) =  0.0; %Es necesario?

        end %for(iElem)
    end
    c_Fext{iSet} = m_Fext(:);
    c_FilFza{iSet} = m_FilFza;

end %for(iSet)

% Ensamble de matriz de fuerzas externas global
Fext = Fext + sparse([c_FilFza{:}],1,cat(1,c_Fext{:}));

    
    

