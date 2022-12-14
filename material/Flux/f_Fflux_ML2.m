function [Fext] = f_Fflux_ML2(Fext,u,e_DatSet,e_VG)
% Determina las fuerzas debidado a las poropresiones en cada paso de tiempo
% y los suma al vector global de fuerzas externas

nSet = e_VG.nSet;
c_Fextxx = cell(nSet,1);
c_FilFzaxx = cell(nSet,1);
c_Fextyy = cell(nSet,1);
c_FilFzayy = cell(nSet,1);
c_Fextzz = cell(nSet,1);
c_FilFzazz = cell(nSet,1);
c_Fextxy = cell(nSet,1);
c_FilFzaxy = cell(nSet,1);
Dtime = e_VG.Dtime;
uxx = u(:,1);
uyy = u(:,2);
uzz = u(:,3);
uxy = u(:,4);
for iSet = 1:nSet
    if e_DatSet(iSet).e_DatMat.conshyp==14||e_DatSet(iSet).e_DatMat.conshyp==15||e_DatSet(iSet).e_DatMat.conshyp==16
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
        uElemSet_exx  = reshape(uxx(dofElemSet),[],nElem);
        uElemSet_eyy  = reshape(uyy(dofElemSet),[],nElem);
        uElemSet_ezz  = reshape(uzz(dofElemSet),[],nElem);
        uElemSet_exy  = reshape(uxy(dofElemSet),[],nElem);
        m_Fextxx = zeros(dofpe,nElem);
        m_Fextyy = zeros(dofpe,nElem);
        m_Fextzz = zeros(dofpe,nElem);
        m_Fextxy = zeros(dofpe,nElem);

        for iElem = 1:nElem
            Hww = zeros(dofpe_p,dofpe_p);
            Lalfa = zeros(dofpe_p,dofpe_p);
            m_Dercae_p = m_DerCa_p(:,:,:,iElem);
            ue_exx = uElemSet_exx(:,iElem);
            ue_eyy = uElemSet_eyy(:,iElem);
            ue_ezz = uElemSet_ezz(:,iElem);
            ue_exy = uElemSet_exy(:,iElem);
            ue_p_exx = ue_exx(pos_p);
            ue_p_eyy = ue_eyy(pos_p);
            ue_p_ezz = ue_ezz(pos_p);
            ue_p_exy = ue_exy(pos_p);

            ue_lambda_exx = ue_exx(pos_lambda);
            ue_lambda_eyy = ue_eyy(pos_lambda);
            ue_lambda_ezz = ue_ezz(pos_lambda);
            ue_lambda_exy = ue_exy(pos_lambda);

            m_pesoPG_p = m_DetJT_p(:,iElem).*wg;
        for iPG = 1:nPG
            N4 = m_FF_p(:,:,iPG);
            DerivN = m_Dercae_p(:,:,iPG);
            Hww = Hww + DerivN'*PermK*DerivN*m_pesoPG_p(iPG); 
            Lalfa = Lalfa + N4'*N4*m_pesoPG_p(iPG);
        end %for(iPG)

        m_Fextxx(pos_p,iElem) =  Dtime*Hww*ue_p_exx;
    %     m_Fextxx(pos_p,iElem) =  Dtime*(Hww*ue_p_exx+Lalfa*ue_lambda_exx);
    %     m_Fextxx(pos_lambda,iElem) =  Dtime*Lalfa*ue_p_exx;
        m_Fextxx(pos_p0,iElem) =  0.0; %Es necesario?
        m_Fextxx(pos_lambda0,iElem) =  0.0; %Es necesario?

        m_Fextyy(pos_p,iElem) =  Dtime*Hww*ue_p_eyy;
    %     m_Fextyy(pos_p,iElem) =  Dtime*(Hww*ue_p_eyy+Lalfa*ue_lambda_eyy);
    %     m_Fextyy(pos_lambda,iElem) =  Dtime*Lalfa*ue_p_eyy;
        m_Fextyy(pos_p0,iElem) =  0.0; %Es necesario?
        m_Fextyy(pos_lambda0,iElem) =  0.0; %Es necesario?

        m_Fextzz(pos_p,iElem) =  Dtime*Hww*ue_p_ezz;
    %     m_Fextzz(pos_p,iElem) =  Dtime*(Hww*ue_p_ezz+Lalfa*ue_lambda_ezz);
    %     m_Fextzz(pos_lambda,iElem) =  Dtime*Lalfa*ue_p_ezz;
        m_Fextzz(pos_p0,iElem) =  0.0; %Es necesario?
        m_Fextzz(pos_lambda0,iElem) =  0.0; %Es necesario?

        m_Fextxy(pos_p,iElem) =  Dtime*Hww*ue_p_exy;
    %     m_Fextxy(pos_p,iElem) =  Dtime*(Hww*ue_p_exy+Lalfa*ue_lambda_exy);
    %     m_Fextxy(pos_lambda,iElem) =  Dtime*Lalfa*ue_p_exy;
        m_Fextxy(pos_p0,iElem) =  0.0; %Es necesario?
        m_Fextxy(pos_lambda0,iElem) =  0.0; %Es necesario?

        end %for(iElem)

        c_Fextxx{iSet} = m_Fextxx(:);
        c_FilFzaxx{iSet} = m_FilFza;

        c_Fextyy{iSet} = m_Fextyy(:);
        c_FilFzayy{iSet} = m_FilFza;

        c_Fextzz{iSet} = m_Fextzz(:);
        c_FilFzazz{iSet} = m_FilFza;

        c_Fextxy{iSet} = m_Fextxy(:);
        c_FilFzaxy{iSet} = m_FilFza;
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
        
        c_Fextxx{iSet} = m_Fext(:);
        c_FilFzaxx{iSet} = m_FilFza;

        c_Fextyy{iSet} = m_Fext(:);
        c_FilFzayy{iSet} = m_FilFza;

        c_Fextzz{iSet} = m_Fext(:);
        c_FilFzazz{iSet} = m_FilFza;

        c_Fextxy{iSet} = m_Fext(:);
        c_FilFzaxy{iSet} = m_FilFza;
    end
end %for(iSet)
% Fext(1:3:end) = -Dtime*Fext(3:3:end); %ESTA BIEN?????????
% Ensamble de matriz de fuerzas externas global
Fext(:,1) = Fext(:,1) + sparse([c_FilFzaxx{:}],1,cat(1,c_Fextxx{:}));
Fext(:,2) = Fext(:,2) + sparse([c_FilFzayy{:}],1,cat(1,c_Fextyy{:}));
Fext(:,3) = Fext(:,3) + sparse([c_FilFzazz{:}],1,cat(1,c_Fextzz{:}));
Fext(:,4) = Fext(:,4) + sparse([c_FilFzaxy{:}],1,cat(1,c_Fextxy{:}));