function m_valHomog = f_HomogArea3d(c_valElem,nVar,volTotal,e_DatSET,e_VG)

   %El espesor o el �rea viene con los pesos de gauss, y por lo tanto no es necesario considerar
   %en la integraci�n (ya viene como peso "volum�trico").

   %%
   %Funci�n que devuelve los valores homogenizados a partir de la integraci�n en toda la superficie
   %de la estructrura.
   %En la celda c_valElem se ingresa las matrices separadas de cada set. Esta matriz debe ser de dimensiones
   %(nVar,npg,nElem) o (nVar*npg,nElem), donde nVar la cantidad de variables o la dimensi�n de la variable que
   %se quiere homogenizar. Esta nVar debe ser la misma para todos los set.
%    nElem = e_VG.nElem;
   nSet = e_VG.nSet;

   m_valHomog = zeros(nVar,1);
   for iSet = 1:nSet
      %nElem = e_DatSET(iSet).nElem;
      wg = e_DatSET(iSet).e_DatElem.wg;
      npg = e_DatSET(iSet).e_DatElem.npg;
      if e_DatSET(iSet).e_DatMat.conshyp==14||e_DatSET(iSet).e_DatMat.conshyp==15||e_DatSET(iSet).e_DatMat.conshyp==16
          c_DetJT = e_DatSET(iSet).m_DetJT_d;
          m_valElemSet = reshape(c_valElem{iSet},nVar,npg,[]);
      elseif e_DatSET(iSet).e_DatMat.conshyp==1
          c_DetJT = e_DatSET(iSet).m_DetJT;
          if isempty(c_valElem{iSet})
              m_valElemSet = 0;
          else
              m_valElemSet = reshape(c_valElem{iSet},nVar,npg,[]);
          end
      end
%       m_valElemSet = reshape(c_valElem{iSet},nVar,npg,[]);
      %
      m_valHomog = m_valHomog+sum(sum(bsxfun(@times,m_valElemSet,...
         reshape(bsxfun(@times,c_DetJT,wg),1,npg,[])),2),3);
   end
   m_valHomog = m_valHomog/volTotal;
      
end
     