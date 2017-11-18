function rotm = calRotationMatrix( v1, v2 )
%calRotationMatrix ������תǰ�������ֵ����ת����
%   rotm ��ת����3*3
%   v1 ��תǰ��������������
%   v2 ��ת���������������

    rotAxis = cross( v1, v2 );
    rotAxis = rotAxis / norm(rotAxis);
    rotAngle = acos( v1 * v2 / (norm(v1) * norm(v2)) );
    % Rodrigues' rotation formula
    K = [          0  -rotAxis(3)   rotAxis(2);
          rotAxis(3)            0  -rotAxis(1);
         -rotAxis(2)   rotAxis(1)           0 ];
    rotm = eye(3) + sin(rotAngle) * K + (1 - cos(rotAngle)) * K^2;
end
