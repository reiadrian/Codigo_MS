function [Fext] = f_Fflux_ML(Fext,u,e_DatSet,e_VG)
% Determina las fuerzas debidado a las poropresiones en cada paso de tiempo
% y los suma al vector global de fuerzas externas

nSet = e_VG.nSet;
c_Fext = cell(nSet,1);
c_FilFza = cell(nSet,1);


for iSet = 1:nSet
    if e_DatSet(iSet).e_DatMat.conshyp==14||e_DatSet(iSet).e_DatMat.conshyp==15||e_DatSet(iSet).e_DatMat.conshyp==16
        Dtime = e_VG.Dtime;
        e_DatMatSet = e_DatSet(iSet).e_DatMat;
        e_DatElemSet = e_DatSet(iSet).e_DatElem;
        m_DerCa_p = e_DatSet(iSet).m_DerCa_p;
        m_FF_p = e_DatSet(iSet).m_FF_p;
        m_DetJT_p = e_DatSet(iSet).m_DetJT_p;
        m_DofElem = e_DatSet(iSet).m_DofElem;
        nElem = e_DatSet(iSet).nElem;

        dofpe = e_DatElemSet.dofpe;
        dofpe_p = e_DatElemSet.dofpe_p;
        wg = e_DatElemSet.wg;
        nPG = e_DatElemSet.npg;
        pos_p = e_DatElemSet.pos_p;
        pos_p0 = e_DatElemSet.pos_p0;
        pos_lambda = e_DatElemSet.pos_lambda;
        pos_lambda0 = e_DatElemSet.pos_lambda0;

        PermK = e_DatMatSet.m_PermK;

        % Grados de libertad y coordenadas de los nodos de los elementos del set
        dofElemSet = m_DofElem(:); 
        m_FilFza = dofElemSet';
        uElemSet  = reshape(u(dofElemSet),[],nElem);

        m_Fext = zeros(dofpe,nElem);
        for iElem = 1:nElem
            Hww = zeros(dofpe_p,dofpe_p);
            Lalfa = zeros(dofpe_p,dofpe_p);
            m_Dercae_p = m_DerCa_p(:,:,:,iElem);
            ue = uElemSet(:,iElem);
            ue_p = ue(pos_p);
            ue_lambda = ue(pos_lambda);
            m_pesoPG_p = m_DetJT_p(:,iElem).*wg;
        for iPG = 1:nPG
            N4 = m_FF_p(:,:,iPG);
            DerivN = m_Dercae_p(:,:,iPG);
            Hww = Hww + DerivN'*PermK*DerivN*m_pesoPG_p(iPG);
            Lalfa = Lalfa + N4'*N4*m_pesoPG_p(iPG);
        end %for(iPG)

    %     m_Fext(pos_p,iElem) =  Dtime*(Hww*ue_p+Lalfa*ue_lambda);
    %     m_Fext(pos_lambda,iElem) =  Dtime*Lalfa*ue_p;
        m_Fext(pos_p,iElem) =  Dtime*Hww*ue_p;
        m_Fext(pos_p0,iElem) =  0.0; %Es necesario?
        m_Fext(pos_lambda0,iElem) =  0.0; %Es necesario?

        end %for(iElem)

        c_Fext{iSet} = m_Fext(:);
        c_FilFza{iSet} = m_FilFza;
    elseif e_DatSet(iSet).e_DatMat.conshyp==1   %Necesito llenar de ceros la     
        e_DatElemSet = e_DatSet(iSet).e_DatElem;
        m_DofElem = e_DatSet(iSet).m_DofElem;
        nElem = e_DatSet(iSet).nElem;
        dofpe = e_DatElemSet.dofpe;
%         wg = e_DatElemSet.wg;              
%         m_DetJT = e_DatSet(iSet).m_DetJT;
%         m_BT = e_DatSet(iSet).m_BT;
%         m_CT = c_CT{iSet};
%         npg = e_DatElemSet.npg;
        % Grados de libertad y coordenadas de los nodos de los elementos del set
        dofElemSet = m_DofElem(:); 
        m_FilFza = dofElemSet';
%         uElemSet  = reshape(u(dofElemSet),[],nElem);
        m_Fext = zeros(dofpe,nElem);
        c_Fext{iSet} = m_Fext(:);
        c_FilFza{iSet} = m_FilFza;
    end
end %for(iSet)

% Ensamble de matriz de fuerzas externas global
Fext = Fext + sparse([c_FilFza{:}],1,cat(1,c_Fext{:}));

    
    

