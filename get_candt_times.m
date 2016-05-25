function [R, r_indx] = get_candt_times(curr_smpl, A)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sample candidate jump times given MJP trajectory
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    prev_t = curr_smpl(1,1);
    prev_s = curr_smpl(2,1);

    R = prev_t;
    for prev_indx = 2:size(curr_smpl,2)
      prev_t_next = curr_smpl(1,prev_indx);
      prev_s_next = curr_smpl(2,prev_indx);

      t_wait = (prev_t_next - prev_t);
      a_s         = A(prev_s);

      n = poissrnd(a_s * t_wait);
      try
        t = rand(n,1).* t_wait + prev_t;
      catch
          n
      end
      R = [R;sort(t); prev_t_next];


      prev_s = prev_s_next;
      prev_t = prev_t_next;
    end;
    assert(issorted(R));
    r_indx = length(R);

