local driver = require"driver"
local image = require"image"
local chronos = require"chronos"

local unpack = table.unpack
local util = require"util"
local blue = require"blue"
local significant = util.significant
local floor, ceil, min, max, sqrt = math.floor, math.ceil, math.min, math.max, math.sqrt

local _M = driver.new()

local B2 = _M.xform(1,-2,1,0,2,-2,0,0,1)
local B3 = {1,-3,3,-1,0,3,-6,3,0,0,3,-3,0,0,0,1}

local function sign(x) return (x<0 and -1) or 1 end
function table.empty (self)
    for _, _ in pairs(self) do  return false end
    return true
end

local function multD(T,n,m,U,r,s) --matrix multiplication T_{mxn} times U_{rxs}
  -- assert(m == r and #T == n*m and #U == r*s, "matrix dimensions dont agree")
  local M = {}
  for i = 1,n do
    for j = 1,s do
      M[s*(i-1)+j] = 0
      for k = 1,m do
        M[s*(i-1)+j] = M[s*(i-1)+j] + T[m*(i-1)+k]*U[s*(k-1)+j]
      end
    end
  end
  return M
end

local function mult(s,u)
  local w = {}
  for i=1,#u do w[i] = s*u[i] end
  return w
end

local function addM(u,v)
  assert(#u == #v, "dot product failed")
  local w = {}
  for i=1,#u do w[i] = u[i]+v[i] end
  return w
end

local function subM(u,v) return addM(u,mult(-1,v)) end

local function dot(u,v)
  assert(#u == #v, "dot product failed")
  return unpack(multD(u,1,#u,v,#v,1))
end

local function norm(v) return dot(v,v) end -- this is the norm squared

local function multP(p,q)
  local n, m, h = #p, #q, {}
  for i=1,(n+m-1) do h[i] = 0 end
  for i=1,n do
    for j=1,m do
      h[i+j-1] = h[i+j-1] + p[i]*q[j]
    end
  end
  return h
end

local function addP(p,q)
  local n, m, h = #p, #q, {}
  for i=1,math.max(m,n) do
    if i <= math.min(m,n) then h[i] = p[i] + q[i]
    elseif n > m then h[i] = p[i]
    else h[i] = q[i]
    end
  end
  return h
end

local function subP(p,q)
  return addP(p,mult(-1,q))
end

local function difP(p)
  local dp = {0}
  for i=1,(#p-1) do dp[i] = i*p[i+1] end
  return dp
end

local function polyToFun(p)
  return function (t)
    local n = (#p-1)
    local v = p[1]
    for i=1,n do v = v + p[i+1]*(t^i) end
    return v
  end
end

local function equal(x,y)
  if x == nil or y == nil then return nil
  else
    return math.abs((x+0.0)-(y+0.0)) < 10^(-8)
  end
end

local function contains(set,x)
  for i,e in pairs(set) do
    if equal(x,e) then return true end
  end
  return false
end

function union(u,v)
  local set = {}
  for i,e in pairs(v) do
    if not contains(set,e) then
      set[#set+1] = e
    end
  end
  for i,e in pairs(u) do
    if not contains(set,e) then
      set[#set+1] = e
    end
  end
  return set
end

local function addRoot(t,s,a,b,roots)
  if t==nil or s==nil then return end
  if sign(s)*t >= sign(s)*s*a and sign(s)*t <= sign(s)*s*b then
    if not contains(roots,t/s) then
      roots[#roots+1] = t/s
    end
  end
  return roots
end

local function newton_aux(f,df,t0,t1,t,dt,y,dy,it)
  if math.abs(dt) < 10^(-12) or it > 50 or y == 0 then
    return t
  else
    if type(dy) ~= "number" or dy == math.huge or ((t-t1)*dy-y)*((t-t0)*dy-y) > 0 or math.abs(2*y) > math.abs(dy*dt) then
      -- bisection
      tm = (t0+t1)*(1/2); dtm = (t1-t0)*(1/2)
      ym = f(tm); dym = df(tm)
      if ym > 0 then
        return newton_aux(f,df,t0,tm,tm,dtm,ym,dym,it+1)
      else
        return newton_aux(f,df,tm,t1,tm,dtm,ym,dym,it+1)
      end
    else
      -- newton step
      dtm = -y*(1/dy); tm = t + dtm
      ym = f(tm); dym = df(tm)
      if ym > 0 then
        return newton_aux(f,df,t0,tm,tm,dtm,ym,dym,it+1)
      else
        return newton_aux(f,df,tm,t1,tm,dtm,ym,dym,it+1)
      end
    end
  end
end

local function safe_newton_raphson(f,df,a,b)
  local y0, y1, x = f(a), f(b), (a+b)*.5
  assert(y0*y1 <= 0, "there is no confirmation of root on the interval")
  if y0 > 0 then
    return newton_aux(f,df,b,a,x,b-a,f(x),df(x),1)
  else
    return newton_aux(f,df,a,b,x,b-a,f(x),df(x),1)
  end
end

local function poly_num_root(p,a,b)
  local Dp = difP(p)
  local funP, fundP = polyToFun(p), polyToFun(Dp)
  local So, S = poly_roots[math.min(#Dp-1,4)](Dp,a,b), {}
  for i,t in pairs(So) do
    if t>a and t<b then S[#S+1] = t end
  end
  local roots = {}
  if funP(a) == 0 then roots[#roots+1] = a end
  if funP(b) == 0 then roots[#roots+1] = b end
  if #S==0 then
    if funP(a)*funP(b) > 0 then return {} end
    return {safe_newton_raphson(funP,fundP,a,b)}
  end
  if funP(a)*funP(S[1]) < 0 then roots[1] = safe_newton_raphson(funP,fundP,a,S[1])
  elseif funP(S[1])==0 then roots[#roots+1] = S[1]
  end
  for i=1,(#S-1) do
    if funP(S[i])*funP(S[i+1]) < 0 then roots[#roots+1] = safe_newton_raphson(funP,fundP,S[i],S[i+1])
    elseif funP(S[i+1]) == 0 then roots[#roots+1] = S[i+1]
    end
  end
  if funP(S[#S])*funP(b) < 0 then roots[#roots+1] = safe_newton_raphson(funP,fundP,S[#S],b)
  elseif funP(S[#S])==0 then roots[#roots+1] = S[#S]
  end
  return roots
end

poly_roots = {
  function (p,a,b)
    local root = -p[1]/p[2]
    if root >= a and root <= b then return {root}
    else return {} end
  end,
  function (p,a,b)
    local a, b = a or -math.huge, b or math.huge
    local n,t1,s1,t2,s2 = _M.quadratic(p[3],p[2],p[1])
    local roots = {}
    addRoot(t1,s1,a,b,roots); addRoot(t2,s2,a,b,roots)
    return roots
  end,
  function (p,a,b)
    local a, b = a or -math.huge, b or math.huge
    local n,t1,s1,t2,s2,t3,s3 = _M.cubic(p[4],p[3],p[2],p[1])
    local roots = {}
    addRoot(t1,s1,a,b,roots); addRoot(t2,s2,a,b,roots); addRoot(t3,s3,a,b,roots)
    return roots
  end,
  function (p,a,b)
    local a, b = a or -math.huge, b or math.huge
    return poly_num_root(p,a,b)
  end
}

local function conicTransf(T,C)
  local A = T.inverse(T)
  return A.transpose(A) * (C * A)
end

local function conicCenter(C)
  return (C[5]*C[3]-C[2]*C[6])/(C[2]^2-C[1]*C[5]), (-C[2]*C[3]+C[1]*C[6])/(C[2]^2-C[1]*C[5])
end

local function applyConic(C,p)
  return dot(p,C * p)
end

local function bezier_to_poly2(b)
  local M = b * B2
  local px, py, pw = {}, {}, {}
  for i=1,3 do
    if #px>0 then table.insert(px,1,M[4-i]) end
    if M[4-i]~=0 and #px==0 then px[1]=M[4-i] end
    if #py>0 then table.insert(py,1,M[7-i]) end
    if M[7-i]~=0 and #py==0 then py[1]=M[7-i] end
    if #pw>0 then table.insert(pw,1,M[10-i]) end
    if M[10-i]~=0 and #pw==0 then pw[1]=M[10-i] end
  end
  return {x=px,y=py,w=pw}
end

local function bezier_to_poly3(b)
  local M = multD(b,2,4,B3,4,4)
  local px, py = {}, {}
  for i=1,4 do
    if #px>0 then table.insert(px,1,M[5-i]) end
    if M[5-i]~=0 and #px==0 then px[1]=M[5-i] end
    if #py>0 then table.insert(py,1,M[9-i]) end
    if M[9-i]~=0 and #py==0 then py[1]=M[9-i] end
  end
  return {x=px,y=py}
end

local function bezier_segment(b,n,t0,t1)
  local levels, newB, t = {b}, {}, 0
  for i=1,(n+1) do
    if i==1 then
      for j=2,n do
        levels[j] = {};
        for k=1,(n+2-j) do
          levels[j][k] = (1-t0)*levels[j-1][k]+t0*levels[j-1][k+1]
          levels[j][n+k+2-j] = (1-t0)*levels[j-1][n+k+3-j]+t0*levels[j-1][n+k+4-j]
          levels[j][2*(n+2-j)+k] = (1-t0)*levels[j-1][2*(n+3-j)+k]+t0*levels[j-1][2*(n+3-j)+k+1]
        end
      end
      newB[1] = (1-t0)*(levels[n][1]) + t0*(levels[n][2])
      newB[n+2] = (1-t0)*(levels[n][3]) + t0*(levels[n][4])
      newB[2*n+3] = (1-t0)*(levels[n][5]) + t0*(levels[n][6])
    else
      if t>0 then
        for j=n+1-t,n do
          for k=1,(n+2-j) do
            levels[j][k] = (1-t1)*levels[j-1][k]+t1*levels[j-1][k+1]
            levels[j][n+k+2-j] = (1-t1)*levels[j-1][n+k+3-j]+t1*levels[j-1][n+k+4-j]
            levels[j][2*(n+2-j)+k] = (1-t1)*levels[j-1][2*(n+3-j)+k]+t1*levels[j-1][2*(n+3-j)+k+1]
          end
        end
      end
      newB[i] = (1-t1)*levels[n][1] + t1*levels[n][2]
      newB[n+i+1] = (1-t1)*levels[n][3] + t1*levels[n][4]
      newB[2*n+i+2] = (1-t1)*levels[n][5] + t1*levels[n][6]
      t = t+1
    end
  end
  return newB
end

local function monotone_partition(p)
  local q = {}
  if p.w == nil then q = multP(difP(p.x),difP(p.y))
  elseif #p.w == 1 then q = multP(difP(p.x),difP(p.y))
  else q = multP(subP(multP(difP(p.x),p.w),multP(p.x,difP(p.w))),subP(multP(difP(p.y),p.w),multP(p.y,difP(p.w))))
  end
  if #q < 2 then return {} end
  local roots = poly_roots[math.min(#q-1,4)](q,0,1)
  table.sort(roots)
  local part, aux = {}, 0
  for i,t in pairs(roots) do
    if t>0 and t<1 and t~=aux then
      part[#part+1] = t
    end
    aux = t
  end
  return part
end

local function cross3(v1,v2,v3)
  local D = {}
  for i = 1,4 do
    local A,s = _M.xform(0,0,0,0,0,0,0,0,0),0
    for j=1,4 do
      if i~=j then
        A[s+1],A[s+2],A[s+3],s = v1[j] or 0,v2[j] or 0,v3[j] or 0,s+3
      end
    end
    D[i] = ((-1)^i)*A.det(A)
  end
  return D
end

local function inflectionPoints(b)
  local D = cross3(b.x,b.y,{1,0,0,0})
  local p = {-D[4], 3*D[3], -3*D[2]}
  local q = {D[3]^2 - D[2]*D[4], - D[2]*D[3], D[2]^2}
  local infs, doublePts = poly_roots[2](p,0,1), poly_roots[2](q,0,1)
  return union(infs,doublePts)
end

local calcWN = {
  function (c,x,y)
    if y > c[2].ymin and y <= c[2].ymax and c[2].l1(x,y) < 0 then
      return c[2].s1
    else return 0 end
  end,
  function (c,x,y)
    if y > c[2].ymin and y <= c[2].ymax then
      if c[2].l1(x,y) <= 0 then
        if c[2].sl > 0 then return c[2].s1
        elseif c[2].R(x,y)*c[2].sR > 0 then return c[2].s1
        end
      else
        if c[2].sl > 0 and c[2].R(x,y)*c[2].sR < 0 then
          return c[2].s1
        end
      end
    end
    return 0
  end,
  function (c,x,y)
    if y > c[2].ymin and y <= c[2].ymax then
      if c[2].l1(x,y) <= 0 then
        if c[2].sl > 0 then return c[2].s1
        else
          if c[2].sT2*c[2].l2(x,y) > 0 and c[2].sT3*c[2].l3(x,y) > 0 then
            if c[2].R(x,y)*c[2].sR > 0 then return c[2].s1 end
          else return c[2].s1
          end
        end
      else
        if c[2].sl > 0 and c[2].sT2*c[2].l2(x,y) >= 0 and c[2].sT3*c[2].l3(x,y) >= 0 then
          if c[2].R(x,y)*c[2].sR < 0 then return c[2].s1 end
        end
      end
    end
    return 0
  end,
  function (c,x,y)
    if applyConic(c[2].conic,_M.vector(x,y)) < 0 then return 1
    else return 0 end
  end
}

local function windingNumber (path,t,x,y)
  local wN = 0
  for i,j in ipairs(t) do
    local c = path[j]
    if c[1]>4 or c[1]<1 then print(c[1]) end
    wN = wN + calcWN[c[1]](c,x,y)
  end
  return wN
end

local function addPoint(box,x,y)
  if table.empty(box) then
    box.xmin, box.xmax, box.ymin, box.ymax = x, x, y, y
  else
    box.xmin, box.xmax, box.ymin, box.ymax = min(box.xmin,x), max(box.xmax,x), min(box.ymin,y), max(box.ymax,y)
  end
  return box
end

local function biP(p,k)
  local l = k%2 + 1
  local n, m, M = #p, 1, {}
  for i=1,n do
    if p[i] then
      M[(m-1)*3 + l] = 0
      M[(m-1)*3 + k] = (i-1)
      M[(m-1)*3 + 3] = p[i]
      m = m+1
    end
  end
  return M
end

local function addbiP(p,q)
  if table.empty(p) then return q
  elseif table.empty(q) then return p
  end
  local h, m, n, i, j, s = {}, #p/3, #q/3, 1, 1, 0
  while  not (i > m and j > n) do
    local px, py, qx, qy = p[(i-1)*3+1] or math.huge, p[(i-1)*3+2] or math.huge, q[(j-1)*3+1] or math.huge, q[(j-1)*3+2] or math.huge
    if (px+py) <= (qx+qy) and (not (px==qx and py==qy)) then
      h[s+1],h[s+2],h[s+3],s,i = px, py, p[i*3],s+3,i+1
    elseif px==qx and py==qy then
      if p[i*3]+q[j*3] then
        h[s+1],h[s+2],h[s+3],s,i,j = px, py,p[i*3]+q[j*3],s+3,i+1,j+1
      end
    else
      h[s+1],h[s+2],h[s+3],s,j = qx, qy, q[j*3],s+3,j+1
    end
  end
  return h
end

local function multbiP(p,q)
  if table.empty(p) then return {0,0,0}
  elseif table.empty(q) then return {0,0,0}
  end
  local h, m, n, order = {}, #p/3, #q/3, {}
  for i=1,m do
    local px, py = p[(i-1)*3+1], p[(i-1)*3+2]
    for j=1,n do
      local qx, qy = q[(j-1)*3+1], q[(j-1)*3+2]
      if h[px+qx]==nil then h[px+qx] = {}; order[#order+1] = px+qx end
      if h[px+qx][py+qy]==nil then h[px+qx][py+qy] = 0 end
      h[px+qx][py+qy] = h[px+qx][py+qy] + q[j*3]*p[i*3]
    end
  end
  table.sort(order)
  local h_, s = {}, 0
  for i,ei in pairs(order) do
    local suborder = {}
    for k in pairs(h[ei]) do suborder[#suborder+1] = k end
    table.sort(suborder)
    for j,ej in pairs(suborder) do
      if h[ei][ej] ~= 0 then
        h_[s+1], h_[s+2], h_[s+3], s = ei, ej, h[ei][ej], s+3
      end
    end
  end
  return h_
end

local function subbiP(p,q)
  return addbiP(p,multbiP({0,0,-1},q))
end

local function det2(M)
  return subbiP(multbiP(M[1],M[4]), multbiP(M[2],M[3]))
end

local function sumbiP(L)
  local total = {0,0,0}
  for i,e in ipairs(L) do
    total = addbiP(total,e)
  end
  return total
end

local function det3(M)
  local a, b, c, d, e, f, g, h, i = unpack(M, 1, 9)
  return sumbiP({multbiP({0,0,-1},multbiP(c,multbiP(e,g))), multbiP(b,multbiP(f,g)), multbiP(c,multbiP(d,h)),
   multbiP({0,0,-1},multbiP(a,multbiP(f,h))), multbiP({0,0,-1},multbiP(b,multbiP(d,i))), multbiP(a,multbiP(e,i))})
end

local function biP_to_fun(p)
  local n = #p
  return function(x,y)
    local v = 0
    for i=1,n,3 do
      local dx, dy, a = p[i], p[i+1], p[i+2]
      v = v + a*(x^dx)*(y^dy)
    end
    return v
  end
end

local implicitEq = {
  function (x1,y1,x2,y2)
    local a, b = y2-y1, x1-x2
    local c, s = -(a*x1 + b*y1), sign(a)
    local p = {0,0,s*c,0,1,s*b,1,0,s*a}
    return biP_to_fun(p), s
  end,
  function(bez)
    local x1,x2,x3 = bez.x[1],bez.x[2] or 0,bez.x[3] or 0
    local y1,y2,y3 = bez.y[1],bez.y[2] or 0,bez.y[3] or 0
    local w1,w2,w3 = 1,0,0
    if bez.w then
      w1,w2,w3 = bez.w[1],bez.w[2] or 0,bez.w[3] or 0
    end
    local f0, f1, f2 = {0,0,-x1,1,0,w1}, {0,0,-x2,1,0,w2}, {0,0,-x3,1,0,w3}
    local g0, g1, g2 = {0,0,-y1,0,1,w1}, {0,0,-y2,0,1,w2}, {0,0,-y3,0,1,w3}
    local bezM = {subbiP(multbiP(f1,g0),multbiP(f0,g1)), subbiP(multbiP(f2,g0),multbiP(f0,g2)),
                          subbiP(multbiP(f2,g0),multbiP(f0,g2)), subbiP(multbiP(f2,g1),multbiP(f1,g2))}
    local p = det2(bezM)
    if p[#p] < 0 then p = multbiP({0,0,-1},p) end
    return biP_to_fun(p)
  end,
  function(bez)
    local x1,x2,x3,x4 = bez.x[1],bez.x[2] or 0,bez.x[3] or 0,bez.x[4] or 0
    local y1,y2,y3,y4 = bez.y[1],bez.y[2] or 0,bez.y[3] or 0,bez.y[4] or 0
    local f0, f1, f2, f3 = {0,0,-x1,1,0,1}, {0,0,-x2}, {0,0,-x3}, {0,0,-x4}
    local g0, g1, g2, g3 = {0,0,-y1,0,1,1}, {0,0,-y2}, {0,0,-y3}, {0,0,-y4}
    local bezM = {subbiP(multbiP(f1,g0),multbiP(f0,g1)),                 subbiP(multbiP(f2,g0),multbiP(f0,g2)),                                subbiP(multbiP(f3,g0),multbiP(f0,g3)),
                  subbiP(multbiP(f2,g0),multbiP(f0,g2)), addbiP(subbiP(multbiP(f3,g0),multbiP(f0,g3)),subbiP(multbiP(f2,g1),multbiP(f1,g2))),  subbiP(multbiP(f3,g1),multbiP(f1,g3)),
                  subbiP(multbiP(f3,g0),multbiP(f0,g3)),                 subbiP(multbiP(f3,g1),multbiP(f1,g3)),                                subbiP(multbiP(f3,g2),multbiP(f2,g3))}
    local p = det3(bezM)
    return biP_to_fun(p)
  end
}

local function intersectLine(l1,l2)
  local a1,b1,a2,b2 = l1[4]-l1[2], l1[1]-l1[3], l2[4]-l2[2], l2[1]-l2[3]
  local c1,c2,s1,s2 = -(a1*l1[1] + b1*l1[2]), -(a2*l2[1] + b2*l2[2]), sign(a1), sign(a2)
  local A = _M.xform(a1,b1,c1,a2,b2,c2,0,0,1)
  if not significant(a1*b2-b1*a2) then return l1[1],l1[2] end
  A = A.inverse(A)
  return A.apply(A,0,0,1)
end

local implicit = {
  function (x1,y1,x2,y2)
    local l1,s1 = implicitEq[1](x1,y1,x2,y2)
    return {xmax = max(x1,x2), ymin = min(y1,y2), ymax = max(y1,y2), l1 = l1, s1 = s1}
  end,
  function(bez, R)
    local x1,y1,x2,y2,x3,y3 = bez[1]/bez[7],bez[4]/bez[7],bez[2]/bez[8],bez[5]/bez[8],bez[3]/bez[9],bez[6]/bez[9]
    local l1, s1 = implicitEq[1](x1,y1,x3,y3)
    local sR,sl = sign(R(x2,y2)), sign(l1(x2,y2))
    return {xmax = max(x1,x3), ymin = min(y1,y3), ymax = max(y1,y3), R = R, l1 = l1, s1 = s1, sR = sR, sl = sl}
  end,
  function(bez, R)
    local x1,y1,x2,y2,x3,y3,x4,y4 = bez[1],bez[5],bez[2],bez[6],bez[3],bez[7],bez[4],bez[8]
    local l1,s1 = implicitEq[1](x1,y1,x4,y4)
    local xT,yT = x1,y1
    if x1==x2 and x1==x3 and y1==y2 and y1==y3 then x2,y2,x3,y3 = x4,y4,x4,y4
    elseif x1==x2 and y1==y2 then x2,y2 = x3,y3
    elseif x4==x3 and y4==y3 then x3,y3 = x2,y2 end
    local l2,s2 = implicitEq[1](x1,y1,x2,y2)
    local l3,s3 = implicitEq[1](x3,y3,x4,y4)
    local xT,yT = intersectLine({x1,y1,x2,y2},{x3,y3,x4,y4})
    if xT==x1 and yT==y1 then R = l1 end
    local sR,sl,sT2,sT3 = sign(R(xT,yT)), sign(l1(xT,yT)), sign(l2(x4,y4)), sign(l3(x1,y1))
    return {xmax = max(x1,x4), ymin = min(y1,y4), ymax = max(y1,y4), R = R, l1 = l1, l2 = l2, l3 = l3, s1 = s1, s2 = s2, s3 = s3, sT2 = sT2, sT3 = sT3, sR = sR, sl = sl}
  end
}

local function transformData(T,data,box,pos,m)
  for i=0,(m-1) do
    data[pos+2*i], data[pos+(2*i+1)] = unpack(T * _M.vector(data[pos+2*i],data[pos+(2*i+1)]))
    addPoint(box,data[pos+2*i],data[pos+2*i+1])
  end
end

-- these are the two functions that you need to modify/implement
local execute = {
  begin_closed_contour = function (shape)
    shape.beg = shape.pos+1
    transformData(shape.xf,shape.data,shape.box,shape.beg,1)
  end,
  begin_open_contour = function (shape)
    shape.beg = shape.pos+1
    transformData(shape.xf,shape.data,shape.box,shape.beg,1)
  end,
  linear_segment = function (shape)
    local data, pos = shape.data, shape.pos
    transformData(shape.xf,data,shape.box,pos+2,1)
    shape.implicitPath[#shape.implicitPath+1] = {1,implicit[1](data[pos],data[pos+1],data[pos+2],data[pos+3])}
  end,
  end_closed_contour = function (shape)
    shape.implicitPath[#shape.implicitPath+1] = {1,implicit[1](shape.data[shape.pos],shape.data[shape.pos+1],shape.data[shape.beg],shape.data[shape.beg+1])}
  end,
  end_open_contour = function (shape)
    if not shape.degenerated then
      shape.implicitPath[#shape.implicitPath+1] = {1,implicit[1](shape.data[shape.pos],shape.data[shape.pos+1],shape.data[shape.beg],shape.data[shape.beg+1])}
    else
      shape.degenerated = false
    end
  end,
  degenerate_segment = function(shape) shape.degenerated = true end,
  quadratic_segment = function (shape)
    local data,pos = shape.data,shape.pos
    transformData(shape.xf,data,shape.box,pos+2,2)
    local b = _M.xform(data[pos],data[pos+2],data[pos+4],data[pos+1],data[pos+3],data[pos+5],1,1,1)
    local partition = monotone_partition(bezier_to_poly2(b))
    table.insert(partition,1,0); partition[#partition+1] = 1
    local R = implicitEq[2](bezier_to_poly2(b));
    for k=1,(#partition-1) do
      local bez = bezier_segment(b,2,partition[k],partition[k+1])
      shape.implicitPath[#shape.implicitPath+1] = {2,implicit[2](bez,R)}
    end
  end,
  cubic_segment = function (shape)
    local data,pos = shape.data,shape.pos
    transformData(shape.xf,data,shape.box,pos+2,3)
    local b = {data[pos],data[pos+2],data[pos+4],data[pos+6],data[pos+1],data[pos+3],data[pos+5],data[pos+7],1,1,1,1}
    local b_poly = bezier_to_poly3(b)
    local partition = union(union(monotone_partition(b_poly),inflectionPoints(b_poly)),{0,1})
    table.sort(partition)
    local d, R = max(#(b_poly.x),#(b_poly.y))-1, {}
    if d>1 then R = implicitEq[d](b_poly)
    else
      R, d = implicitEq[1](b[1],b[5],b[4],b[8]), 1
    end
    for k=1,(#partition-1) do
      local bez = bezier_segment(b,3,partition[k],partition[k+1])
      shape.implicitPath[#shape.implicitPath+1] = {d,implicit[3](bez,R)}
    end
  end,
  rational_quadratic_segment = function (shape)
    local data,pos = shape.data,shape.pos
    data[pos+2], data[pos+3], data[pos+4] = unpack(multD(shape.xf,3,3,{data[pos+2], data[pos+3], data[pos+4]},3,1))
    data[pos+5], data[pos+6] = unpack(multD(shape.xf,3,3,{data[pos+5], data[pos+6], 1},3,1))
    addPoint(shape.box,data[pos+5], data[pos+6]); addPoint(shape.box,data[pos+2]*(1/data[pos+4]), data[pos+3]*(1/data[pos+4]))
    local b = _M.xform(data[pos],data[pos+2],data[pos+5],data[pos+1],data[pos+3],data[pos+6],1,data[pos+4],1)
    local partition = monotone_partition(bezier_to_poly2(b))
    table.insert(partition,1,0); partition[#partition+1] = 1
    local R = implicitEq[2](bezier_to_poly2(b));
    for k=1,(#partition-1) do
      local bez = bezier_segment(b,2,partition[k],partition[k+1])
      shape.implicitPath[#shape.implicitPath+1] = {2,implicit[2](bez,R)}
    end
  end
}

local spread = {
  pad = function (t)
    return min(1,max(0,t))
  end,
  ["repeat"] = function (t)
    return t-floor(t)
  end,
  reflect = function (t)
    return 2*math.abs(t*.5 - floor(t*.5 + .5))
  end,
  transparent = function (t) return t end
}

local function interpolate(t1,t2,t,c1,c2)
  return mult(1/(t2-t1),addM(mult(t2-t,c1),mult(t-t1,c2)))
end

local function ramp(r,t)
  for i=1,(#r-3),2 do
    if t >= r[i] and t < r[i+2] then
      return interpolate(r[i],r[i+2],t,r[i+1],r[i+3])
    end
  end
  if r.spread == "transparent" then return {0,0,0,0}
  elseif t < r[1] then return r[2]
  elseif t >= r[#r-1] then return r[#r]
  end
  return {0,0,0,0}
end

local function ramp_texture(img,p)
  local w, h = img.width, img.height
  local x, y = p[1]*w + 1/2, p[2]*h + 1/2
  local x1,y1,x2,y2,x3,y3,x4,y4 = floor(x), ceil(y), ceil(x), ceil(y), floor(x), floor(y), ceil(x), floor(y)
  local a1,a2,a3,a4 = (x4-x)*(y-y4), (x-x3)*(y-y3), (x2-x)*(y2-y), (x-x1)*(y1-y)
  local c1,c2,c3,c4 = {img:get(max(x1,1),min(y1,h))}, {img:get(min(x2,w),min(y2,h))}, {img:get(max(x3,1),max(y3,1))}, {img:get(min(x4,w),max(y4,1))}
  return addM(addM(mult(a1,c1),mult(a2,c2)),addM(mult(a3,c3),mult(a4,c4)))
end


local prepare = {
  start = function(scene, viewport)
    scene.width, scene.height = viewport[3]-viewport[1], viewport[4]-viewport[2]
    local n = floor(max(scene.width,scene.height)*(1/20))
    scene.grid = {}
    for i = 1,n do scene.grid[i] = {} for j = 1,n do scene.grid[i][j] = {} end end --creating a grid
  end,
  triangle = function (shape)
    local x1,x2,x3,y1,y2,y3 = unpack(shape.xf * _M.xform(shape.x1,shape.x2,shape.x3,shape.y1,shape.y2,shape.y3,1,1,1))
    shape.box = {xmin = min(x1,x2,x3), xmax = max(x1,x2,x3), ymin = min(y1,y2,y3), ymax = max(y1,y2,y3)}
    shape.implicitPath = {{1,implicit[1](x1,y1,x2,y2)},{1,implicit[1](x2,y2,x3,y3)},{1,implicit[1](x3,y3,x1,y1)}}
  end,
  polygon = function (shape)
    local data = shape.data
    local n = #data
    shape.implicitPath, shape.box = {}, {}
    transformData(shape.xf, data, shape.box, 1, n*.5)
    for i=1,(n-3),2 do
      shape.implicitPath[#shape.implicitPath+1] =  {1,implicit[1](data[i],data[i+1],data[i+2],data[i+3])}
    end
    shape.implicitPath[#shape.implicitPath+1] = {1,implicit[1](data[n-1],data[n],data[1],data[2])}
  end,
  circle = function (shape)
    assert(shape.r > 0, "invalid value for radius")
    local C = conicTransf(shape.xf, _M.xform(1,0,-shape.cx,0,1,-shape.cy,-shape.cx,-shape.cy,shape.cx^2+shape.cy^2-shape.r^2))
    local cx, cy= conicCenter(C)
    local C_ = conicTransf(_M.affine(1,0,-cx,0,1,-cy),C)
    local ry = math.sqrt(C_[1]*C_[9]*(1/(C_[2]^2 - C_[1]*C_[5]))) --calculating conic bounding box
    local rx = math.sqrt(C_[5]*C_[9]*(1/(C_[2]^2 - C_[1]*C_[5])))
    shape.box = {xmin = cx-rx, xmax = cx+rx, ymin = cy-ry, ymax = cy+ry}
    shape.implicitPath = {{4,{xmax = shape.box.xmax, ymin = shape.box.ymin,ymax = shape.box.ymax,conic = C}}}
  end,
  path = function (shape)
    shape.implicitPath, shape.box = {}, {}
    shape.beg, shape.degenerated = 1, false
    for j, inst in pairs(shape.instructions) do
      shape.pos = shape.offsets[j]
      execute[inst](shape)
    end
  end,
  grid = function (scene,viewport,shape,i)
    local x0,y0,xn,yn = unpack(viewport, 1, 4)
    local n = #scene.grid
    local function indexBound(x,t0,tn)
      local w = n*(x-t0)*(1/(tn-t0))
      if w <= 0 then return 1
      elseif w >= n then return n
      else return ceil(w) end
    end
    local i1 = indexBound(shape.box.xmin,x0,xn) -- Change i1 for each curve
    for j = 1,#shape.implicitPath do
      local i2, j1, j2 = indexBound(shape.implicitPath[j][2].xmax,x0,xn), indexBound(shape.implicitPath[j][2].ymin,y0,yn), indexBound(shape.implicitPath[j][2].ymax,y0,yn)
      for k = i1,i2 do
        for l = j1,j2 do
          scene.grid[k][l][i] = scene.grid[k][l][i] or {}
          scene.grid[k][l][i][#scene.grid[k][l][i]+1] = j
        end
      end
    end
  end,
  conclude = function(scene, viewport)
    local n = #scene.grid
    for i = 1,n do for j = 1,n do
      local order = {}
      local size = 0
      for k,e in pairs(scene.grid[i][j]) do
        order[#order+1] = k
        size = size + #e
      end
      table.sort(order)
      scene.grid[i][j].order = order
      scene.grid[i][j].size = size
    end end --sorting grid
    scene.vxmin, scene.vymin  = viewport[1], viewport[2]
  end,
  solid = function (paint)
    paint.data[4] = paint.data[4]*paint.opacity
    paint.color = function (self,x,y)
      return self.data
    end
  end,
  lineargradient = function (paint)
    local data = paint.data
    for i=1,#data.ramp-1,2 do data.ramp[i+1][4] = data.ramp[i+1][4]*paint.opacity end
    local wrapping = spread[data.ramp.spread]
    local function lmap(p)
      local v = data.p2 - data.p1
      return dot(p-data.p1,v)/norm(v) -- a lot of divisions :/
    end
    paint.color = function (self,x,y)
      return ramp(data.ramp,wrapping(lmap(self.xf*_M.vector(x,y))))
    end
  end,
  radialgradient = function (paint)
    local data = paint.data
    for i=1,#data.ramp-1,2 do data.ramp[i+1][4] = data.ramp[i+1][4]*paint.opacity end
    local r,c,f = data.radius, data.center, data.focus
    assert(r > 0, "wrong parameters")
    local wrapping = spread[data.ramp.spread]
    local v = f-c
    local d_to_v = norm(v)
    if d_to_v > r^2 then
      r = c + (r/math.sqrt(d_to_v))*v
    end
    local function rmap(p)
      local u = p-f
      local d1 = norm(u)
      if d1 == 0 then return 0 end
      local p = {d_to_v-r^2,2*dot(u,v),d1}
      local t = unpack(poly_roots[2](p,0,math.huge))
      if t == 0 or t == nil then return 1
      else return 1/t end
    end
    paint.color = function (self,x,y)
      return ramp(data.ramp,wrapping(rmap(self.xf*_M.vector(x,y))))
    end
  end,
  texture = function (paint)
    local img = paint.data.image
    local wrapping = function(x,y) return _M.vector(spread.pad(x), spread.pad(y)) end
    paint.color = function (self,x,y)
      return ramp_texture(img,wrapping(unpack(self.xf*_M.vector(x,y))))
    end
  end
}

-- prepare scene for sampling and return modified scene
local function preparescene(scene, viewport)
  prepare.start(scene, viewport)
  for i, element in ipairs(scene.elements) do
    element.shape.xf = scene.xf * element.shape.xf
    prepare[element.shape.type](element.shape)
    prepare["grid"](scene, viewport, element.shape, i)
    element.paint.xf = element.paint.xf:inverse() * scene.xf:inverse()
    prepare[element.paint.type](element.paint)
  end
  -- print(totalTime)
  prepare.conclude(scene, viewport)
  return scene
end

local function inside(box,x,y)
  return x >= box.xmin and x <= box.xmax and y >= box.ymin and y <= box.ymax
end

local function blend(f,b)
  return {f[4]*f[1]+(1-f[4])*b[4]*b[1],f[4]*f[2]+(1-f[4])*b[4]*b[2],f[4]*f[3]+(1-f[4])*b[4]*b[3],f[4]+(1-f[4])*b[4]}
end

local function ungamma(c)
  return {c[1]^2,c[2]^2,c[3]^2,c[4]}
end

local function gamma(c)
  local c1,c2,c3,c4 = unpack(c)
  return {sqrt(c1),sqrt(c2),sqrt(c3),c4}
end
-- sample scene at x,y and return r,g,b,a
local function sample(scene, x, y)
  local sum = {0,0,0,0}
  local k,l = ceil((#scene.grid)*(x-scene.vxmin)*(1/scene.width)), ceil((#scene.grid)*(y-scene.vymin)*(1/scene.height))
  local m, n_samples = scene.grid[k][l].size, 1
  if m >= 5 and m < 15 then n_samples = 8
  elseif m >= 15 and m < 25 then n_samples = 16
  elseif m >= 35 and m < 55 then n_samples = 32
  elseif m >= 55 then n_samples = 64
  end
  for j = 1,n_samples do
    local scolor = {1,1,1,1}
    local sx, sy = x + blue[n_samples][j], y + blue[n_samples][j+1]
    for ind,i in ipairs(scene.grid[k][l].order) do
      local t = scene.grid[k][l][i]
      local shape, paint = scene.elements[i].shape, scene.elements[i].paint
      if inside(shape.box,sx,sy) then
        local wN = windingNumber(shape.implicitPath,t,sx,sy)
        if scene.elements[i].type == "fill" and wN ~= 0 then
          scolor = blend(paint:color(sx,sy),scolor)
        elseif scene.elements[i].type == "eofill" and wN % 2 == 1 then
          scolor = blend(paint:color(sx,sy),scolor)
        end
      end
    end
    sum = addM(sum,ungamma(scolor))
  end
  color = gamma(mult(1/n_samples,sum))
  return unpack(color,1,4)
end

-- verifies that there is nothing unsupported in the scene
local function checkscene(scene)
    for i, element in ipairs(scene.elements) do
        assert(element.type == "fill" or
               element.type == "eofill", "unsupported element")
        assert(element.shape.type == "circle" or
               element.shape.type == "triangle" or
               element.shape.type == "path" or
               element.shape.type == "polygon", "unsuported primitive")
        assert(element.paint.type == "solid" or
               element.paint.type == "lineargradient" or
               element.paint.type == "texture" or
               element.paint.type == "radialgradient", "unsupported paint")
    end
end

-- output formatted string to stderr
local function stderr(...)
    io.stderr:write(string.format(...))
end

function _M.render(scene, viewport, file)
local time = chronos.chronos()
    -- make sure scene does not contain any unsuported content
    checkscene(scene)
    -- transform and prepare scene for rendering
    scene = preparescene(scene, viewport)
    -- get viewport
    local vxmin, vymin, vxmax, vymax = unpack(viewport, 1, 4)
stderr("preprocess in %.3fs\n", time:elapsed())
time:reset()
    -- get image width and height from viewport
    local width, height = vxmax-vxmin, vymax-vymin
    -- allocate output image
    local img = image.image(width, height)
    -- render
    for i = 1, height do
stderr("\r%5g%%", floor(1000*i/height)/10)
        local y = vymin+i-1.+.5
        for j = 1, width do
            local x = vxmin+j-1.+.5
            img:set(j, i, sample(scene, x, y))
        end
    end
stderr("\n")
stderr("rendering in %.3fs\n", time:elapsed())
time:reset()
    -- store output image
    image.png.store8(file, img)
stderr("saved in %.3fs\n", time:elapsed())
end

return _M
