local matrix ={}

local mat_obj = {
    data = {},
    __type = "matrix"
}
local mat_meta = {}

local esp = 0.0001

local function deep_copy(object)
    local lookup_table = {}
    local function _copy(obj)
        if type(obj) ~= "table" then
            return obj
        elseif lookup_table[obj] then
            return lookup_table[obj]
        end
        local new_table = {}
        lookup_table[obj] = new_table
        for key, value in pairs(obj) do
            new_table[_copy(key)] = _copy(value)
        end
        return setmetatable(new_table, getmetatable(obj))
    end
    return _copy(object)
end

local function is_mat(mat)
    if type(mat) ~= "table" then
        return false
    end
    if type(mat.type) ~= "function" then
        return false
    end
    if mat:type() ~= "matrix" then
        return false
    end
    return true
end

---------------- init the matrix --------------------
local function mat_init()
    return deep_copy(mat_obj)
end

function matrix.zeros(row, column)
    local mat = mat_init()
    for i = 1, row do
        mat.data[i] = {}
        for j = 1, column do
            mat.data[i][j] = 0
        end
    end
    setmetatable(mat, mat_meta)
    return mat
end

function matrix.I(m)
    local mat = matrix.zeros(m, m)
    for i = 1, m do
        mat.data[i][i] = 1
    end
    return mat
end

function matrix.new(row, ...)
    local rows = {row, ...}
    local column = #row
    local ret = matrix.zeros(#rows, column)
    for i = 1, #rows, 1 do
        local r = rows[i]
        if column ~= #r then
            print("[ERROR] At least one row has different column than others.")
            return nil
        end
        for j = 1, column, 1 do
            local c = r[j]
            if type(c) ~= "number" then
                print("[ERROR] At least one element is NOT type number.")
                return nil
            else
                ret.data[i][j] = c
            end
        end
    end
    return ret
end
-----------------------------------------------------

------------- base function of the matrix -----------
function mat_obj:type()
    return tostring(self.__type)
end

function mat_obj:size()
    local row, column = #self.data, #self.data[1]
    return row, column
end

function mat_meta:__tostring()
    local row, column = self:size()
    local ret = "["
    for i = 1, row do
        if i ~= 1 then
            ret = ret .. " "
        end
        ret = ret .. "["
        for j = 1, column do
            ret = ret .. tostring(self.data[i][j])
            if j ~= column then
                ret = ret .. " "
            end
        end
        ret = ret .. "]"
        if i ~= row and column ~= 1 then
            ret = ret .. "\n"
        end
    end
    ret = ret .. "]"
    return ret
end

function mat_obj:print()
    print(tostring(self))
end

function mat_obj:get(row, column)
    return self.data[row][column]
end

function mat_meta:__add(mat)
    if is_mat(mat) == false then
        print("[ERROR] matrix can add / sub matrix ONLY.")
        return nil
    end
    local self_row, self_column = self:size()
    local mat_row, mat_column = mat:size()
    if self_row ~= mat_row or self_column ~= mat_column  then
        print(string.format(
                "[ERROR] matrix <%dx%d> cannot add / sub with <%dx%d>",
                self_row, self_column, mat_row, mat_column
            )
        )
        return nil
    end
    local ret = matrix.zeros(self_row, self_column)
    for i = 1, self_row, 1 do
        for j = 1, self_column, 1 do
            ret.data[i][j] = self.data[i][j] + mat.data[i][j]
        end
    end
    return ret
end

function mat_meta:__mul(mat_or_num)
    if is_mat(mat_or_num) == true then
        local self_row, self_column = self:size()
        local mat_row, mat_column = mat_or_num:size()
        if self_column ~= mat_row  then
            print(string.format(
                    "[ERROR] matrix <%dx%d> cannot multiply with <%dx%d>",
                    self_row, self_column, mat_row, mat_column
                )
            )
            return nil
        end
        local ret = matrix.zeros(self_row, mat_column)
        for i = 1, self_row, 1 do
            for j = 1, mat_column, 1 do
                for k = 1, self_column, 1 do
                    ret.data[i][j] = ret.data[i][j] + self.data[i][k] * mat_or_num.data[k][j]
                end
            end
        end
        return ret
    elseif type(mat_or_num) == "number" then
        local ret = deep_copy(self)
        local self_row, self_column = self:size()
        for i = 1, self_row, 1 do
            for j = 1, self_column, 1 do
                ret.data[i][j] = ret.data[i][j] * mat_or_num
            end
        end
        return ret
    else
        print("[ERROR] matrix can multiply matrix or number ONLY.")
        return nil
    end
end

function mat_meta:__sub(mat)
    return self + mat * -1
end

function mat_meta:__eq(mat)
    if is_mat(mat) == false then
        print("[ERROR] matrix can compare with matrix ONLY.")
        return nil
    end
    local self_row, self_column = self:size()
    local mat_row, mat_column = mat:size()
    if self_row ~= mat_row or self_column ~= mat_column  then
        return false
    else
        for i= 1, self_row, 1 do
            for j = 1, self_column, 1 do
                if math.abs(self.data[i][j] - mat.data[i][j]) > esp then
                    return false
                end
            end
        end
        return true
    end
end

function mat_obj:T()
    local self_row, self_column = self:size()
    local ret = matrix.zeros(self_column, self_row)
    for i = 1, self_column, 1 do
        for j = 1, self_row, 1 do
            ret.data[i][j] = self.data[j][i]
        end
    end
    return ret
end

-- ref: https://rosettacode.org/wiki/LU_decomposition#JavaScript
function mat_obj:pivot()
    local self_row, self_column = self:size()
    if self_row ~= self_column then
        print(string.format(
            "[ERROR] matrix <%dx%d> does NOT have pivotized matrix",
            self_row, self_column)
        )
        return nil
    end
    local ret = matrix.I(self_row)
    for i = 1, self_column, 1 do
        local maxm = self.data[i][i]
        local curr_row = i
        for j = i, self_row, 1 do
            if self.data[j][i] > maxm then
                maxm = m[j][i]
                row = j
            end
        end
        if i ~= curr_row then   -- switch the lines
            local tmp = ret.data[i]
            ret.data[i] = ret.data[row]
            ret.data[row] = tmp
        end
    end
    return ret
end

function mat_obj:LUPdecompose()
    -- PLU decompostion solution
    -- ref: https://en.wikipedia.org/wiki/LU_decomposition
    local self_row, self_column = self:size()
    if self_row ~= self_column then
        print(string.format(
            "[ERROR] matrix <%dx%d> cannot be L-U decomposed.",
            self_row, self_column)
        )
        return nil
    end
    local A = deep_copy(self)
    local L = matrix.I(self_row)
    local U = matrix.zeros(self_row, self_column)
    local P = matrix.I(self_row)
    local det_P = 1
    for i = 1, self_row, 1 do
        -- finding the pivoting row
        local maxNum = 0.0
        local rowID_2_find_maxNum = i
        for k = i, self_row, 1 do -- find the below rows
            local absNum = math.abs(A.data[k][i])
            --print("absNum", absNum)
            if absNum > maxNum then
                maxNum = absNum
                rowID_2_find_maxNum = k
            end
        end
        if maxNum < esp then
            local failed = false
            if i ~= self_row then
                failed = true
            else
                for c = 1, self_column, 1 do
                    if math.abs(A.data[i][c]) > esp then
                        failed = false
                        break
                    else
                        failed = true
                    end
                end
            end
            if failed then
                print("[ERROR] this matrix is degenerate.")
                return nil
            end
        end
        -- switching the row
        if rowID_2_find_maxNum ~= i then
            -- pivoting P
            local tmp_row = P.data[i]
            P.data[i] = P.data[rowID_2_find_maxNum]
            P.data[rowID_2_find_maxNum] = tmp_row
            -- pivoting matrix A
            tmp_row = A.data[i]
            A.data[i] = A.data[rowID_2_find_maxNum]
            A.data[rowID_2_find_maxNum] = tmp_row
            -- calc det_P
            det_P = det_P * -1
        end
        for k = i, self_column, 1 do
            local sum = 0
            for j = 1, i, 1 do
                sum = sum + (L.data[i][j] * U.data[j][k])
            end
            U.data[i][k] = A.data[i][k] - sum
        end
        for k = i, self_column, 1 do
            if i == k then
                L.data[i][i] = 1;
            else
                local sum = 0
                for j = 1, i, 1 do
                    sum = sum + L.data[k][j] * U.data[j][i]
                end
                L.data[k][i] = (A.data[k][i] - sum) / U.data[i][i]
            end
        end
    end
    return L, U, P, det_P
end

--function mat_obj:determinant()
    --local ret = 0
    --local row, column = self:size()
    --if row ~= column then
        --print(string.format("[ERROR] matrix <%sx%s> does NOT have determinant", row, column))
        --return nil
    --end
--end

local function lower_triangle_inverse(l)
    local n, _ = l:size()
    local ret = matrix.zeros(n, n)
    for i = 1, n, 1 do
        ret.data[i][i] = 1 / l.data[i][i]
        for j = 1, i - 1, 1 do
            local s = 0
            for k = j, i - 1, 1 do
                s = s + l.data[i][k] * ret.data[k][j]
            end
            ret.data[i][j] = -s * ret.data[i][i]
        end
    end
    return ret
end

-- since PA = P * L * U, thus: A^-1 = U^-1 * L^-1 * P^-1, and P^T = P^-1
-- matrix U^-1 & L^-1 can be caculated by Elementary row operations
function mat_obj:inverse()
    local L, U, P, _ = self:LUPdecompose()
    local self_row, self_column = self:size()
    if self_row ~= self_column then
        print(string.format(
            "[ERROR] matrix <%dx%d> does NOT have inversed matrix",
            self_row, self_column)
        )
        return nil
    end
    local l_inverse = lower_triangle_inverse(L)
    local u_inverse = lower_triangle_inverse(U:T()):T()
    return u_inverse * l_inverse * P:T()
end
-----------------------------------------------------

return matrix
