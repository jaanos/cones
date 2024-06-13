# used in a subset of the polyhedral package
# to use in a standard install of polyhedral, change
# polyhedral-mini -> polyhedral
# Janoš Vidali 2010
FileSCDD:=Filename(DirectoriesPackagePrograms("polyhedral-mini"),"scdd_gmp");

OnMultiPairs := function(x, g)
    if x[1] = x[2] then
        return OnTuples(x, g);
    else
        return OnSets(x, g);
    fi;
end;

OnIndices := function(x, g)
    local h;
    h := g^-1;
    return List([1..Length(x)], i -> x[OnPoints(i, h)]);
end;

MakeIndices := function(n, ord, diag)
    local L, i, j, k, f;
    L := [];
    for i in [1..n] do
        if ord then
            k := n;
        else
            k := i;
        fi;
        for j in [1..k] do
            if j <> i or diag then
                Add(L, [j, i]);
            fi;
        od;
    od;
    return L;
end;

IndexSymmGroup := function (n, L, ord)
    local G, f;
    G := SymmetricGroup(n);
    if ord then
        f := OnTuples;
    else
        f := OnMultiPairs;
    fi;
    return Action(G, L, f);
end;

ExpandRepr := function (repr, G)
    local L, orb, len, pos, r;
    Print("Expanding representation with ", String(Length(repr)), " orbits...\n");
    L := [];
    len := [];
    pos := [];
    for r in repr do
        orb := Orbit(G, r, OnIndices);
        Add(pos, Length(L)+1);
        Append(L, orb);
        Add(len, Length(orb));
    od;
    Print("Expanded representation has ", String(Length(L)), " generators.\n");
    return rec(L := L, pos := pos, len := len);
end;

CollapseRepr := function (repr, G, inc)
    local L, o, r, i, red;
    Print("Collapsing representation with ", String(Length(repr)), " generators...\n");
    L := [];
    red := [];
    i := 0;
    for r in repr do
        i := i+1;
        if r in L then
            continue;
        fi;
        o := Orbit(G, r, OnIndices);
        Add(red, rec(repr := r, inc := inc[i], len := Length(o)));
        Append(L, o);
    od;
    Print("Collapsed representation has ", String(Length(red)), " orbits.\n");
    return red;
end;

SortByIncidence := function (x, y)
    if x.inc <> y.inc then
        return x.inc > y.inc;
    elif x.len <> y.len then
        return x.len > y.len;
    else
        return x.repr < y.repr;
    fi;
end;

SortByLen := function (x, y)
    if x.len <> y.len then
        return x.len > y.len;
    else
        return x.p < y.p;
    fi;
end;

ReadNormalized := function (f)
    local l;
    l := ReadLine(f);
    if l = fail then
        return fail;
    fi;
    NormalizeWhitespace(l);
    return l;
end;

ReadRepr := function (filename)
    local a, o, m, n, f, l, t, M, N, L, stage, ord, diag, done, vrep, hrep;
    Print("Reading file ", filename, " ...\n");
    ord := false;
    diag := false;
    vrep := false;
    hrep := false;
    done := true;
    stage := 0;
    M := [];
    N := [];
    L := [];
    f := InputTextFile(filename);
    while true do
        l := ReadNormalized(f);
        if l = fail then
            break;
        fi;
        if IsEmptyString(l) then
            continue;
        fi;
        if stage = 0 then
            if l = "V-representation" then
                vrep := true;
            elif l = "H-representation" then
                hrep := true;
            elif l = "ordered" then
                ord := true;
            elif l = "unordered" then
                ord := false;
            elif l = "diagonal" then
                diag := true;
            elif l = "no diagonal" then
                diag := false;
            elif l = "begin" then
                if not vrep and not hrep then
                    Print("ERROR: none of V- and H- representation specified!");
                    done := false;
                    break;
                fi;
                stage := 1;
                continue;
            fi;
            if vrep and hrep then
                Print("ERROR: both V- and H- representation specified!");
                done := false;
                break;
            fi;
        elif stage = 1 then
            a := SplitString(l, " ");
            o := Int(a[1]);
            m := Int(a[2]);
            t := a[3];
            n := m*(m-1)/2;
            if ord then
                n := 2*n;
            fi;
            if diag then
                n := n + m;
            fi;
            stage := 2;
        elif stage = 2 then
            L := List(SplitString(l, " "){[1..n]},
                function (x)
                    local z;
                    z := List(SplitString(x, ","), Int);
                    if not ord then
                        Sort(z);
                    fi;
                    return z;
                end);
            stage := 3;
        elif stage = 3 then
            if l = "end" then
                stage := 4;
                continue;
            fi;
            Add(M, SplitString(l, " "){[1..n]});
        elif stage = 4 then
            if l = "nullspace" then
                stage := 5;
            fi;
        elif stage = 5 then
            stage := 6;
        elif stage = 6 then
            if l = "end" then
                break;
            fi;
            Add(N, SplitString(l, " "){[1..n]});
        fi;
    od;
    CloseStream(f);
    if not done then
        return fail;
    fi;
    return rec(m := m,  n:= n, L := L, M := M, N := N, t := t, ord := ord, diag := diag, hrep := hrep);
end;

WriteCDD := function(basename, n, R, t, hrep)
    local cddinput, cddoutput, cddinc, f, l;
    
    if hrep then
        cddinput := Concatenation(basename, ".ine");
        cddoutput := Concatenation(basename, ".ext");
        cddinc := Concatenation(basename, ".ecd");
    else
        cddinput := Concatenation(basename, ".ext");
        cddoutput := Concatenation(basename, ".ine");
        cddinc := Concatenation(basename, ".icd");
    fi;
    Print("Writing file ", cddinput, " ...\n");
    f := OutputTextFile(cddinput, false);
    
    if hrep then
        WriteLine(f, "H-representation");
    else
        WriteLine(f, "V-representation");
    fi;
    WriteLine(f, "begin");
    WriteLine(f, Concatenation(" ", String(Length(R)), " ", String(n+1), " ", t));
    
    for l in R do
        WriteLine(f, Concatenation(" 0 ", JoinStringsWithSeparator(l, " ")));
    od;
    
    WriteLine(f, "end");
    CloseStream(f);
    return rec(cddinput := cddinput, cddoutput := cddoutput, cddinc := cddinc);
end;

ReadCDD := function (filename)
    local l, f, a, n, t, R, stage;
    Print("Reading file ", filename, " ...\n");
    R := [];
    stage := 0;
    f := InputTextFile(filename);
    while true do
        l := ReadNormalized(f);
        if l = fail or l = "end" then
            break;
        fi;
        if IsEmptyString(l) then
            continue;
        fi;
        if stage = 0 then
            if l = "begin" then
                stage := 1;
            fi;
        elif stage = 1 then
            a := SplitString(l, " ");
            n := Int(a[2]);
            t := a[3];
            stage := 2;
        elif stage = 2 then
            Add(R, SplitString(l, " "){[2..n]});
        fi;
    od;
    CloseStream(f);
    return rec(R := R, t := t);
end;

ReadCDDInc := function (filename)
    local l, f, inc, stage;
    Print("Reading file ", filename, " ...\n");
    inc := [];
    stage := 0;
    f := InputTextFile(filename);
    while true do
        l := ReadNormalized(f);
        if l = fail or l = "end" then
            break;
        fi;
        if IsEmptyString(l) then
            continue;
        fi;
        if stage = 0 then
            if l = "begin" then
                stage := 1;
            fi;
        elif stage = 1 then
            stage := 2;
        elif stage = 2 then
            Add(inc, AbsInt(Int(SplitString(l, " ")[2])));
        fi;
    od;
    CloseStream(f);
    return inc;
end;

WriteRepr := function(basename, m, L, M, nEXT, t, ord, diag, hrep)
    local filename, f, l, fields;
    
    if hrep then
        filename := Concatenation(basename, ".oin");
    else
        filename := Concatenation(basename, ".oex");
    fi;
    Print("Writing file ", filename, " ...\n");
    f := OutputTextFile(filename, false);
    
    if hrep then
        WriteLine(f, "oin_file: Orbits");
        WriteLine(f, Concatenation("name ", filename));
        WriteLine(f, "H-representation");
    else
        WriteLine(f, "oex_file: Orbits");
        WriteLine(f, Concatenation("name ", filename));
        WriteLine(f, "V-representation");
    fi;
    
    if ord then
        WriteLine(f, "ordered");
    else
        WriteLine(f, "unordered");
    fi;
    
    if diag then
        WriteLine(f, "diagonal");
    else
        WriteLine(f, "no diagonal");
    fi;
    
    for l in M do
        fields := RecFields(l);
        if not ("inc" in fields) then
            l.inc := "";
        fi;
        if not ("len" in fields) then
            l.len := "";
        fi;
    od;
    Sort(M, SortByIncidence);
    
    WriteLine(f, Concatenation("total ", String(Sum(List(M, x -> x.len)))));
    WriteLine(f, "begin");
    WriteLine(f, Concatenation("\t", String(Length(M)), " ", String(m), " ", t));
    WriteLine(f, Concatenation("\t", JoinStringsWithSeparator(List(L,
        x -> JoinStringsWithSeparator(List(x, String), ",")), "\t"), "\tinc\ts"));
    
    for l in M do
        WriteLine(f, Concatenation("\t", JoinStringsWithSeparator(l.repr, "\t"),
            "\t", String(l.inc), "\t", String(l.len)));
    od;
    
    WriteLine(f, "end");
    
    if Length(nEXT) > 0 then
        WriteLine(f, "nullspace");
        WriteLine(f, Concatenation("\t", String(Length(nEXT))));
        for l in nEXT do
            WriteLine(f, Concatenation("\t", JoinStringsWithSeparator(l, "\t")));
        od;
        WriteLine(f, "end");
    fi;
    
    CloseStream(f);
end;

ConvertRepr := function (filename)
    local basename, repr, G, R, write, read, inc, M;
    basename := SplitString(filename, ".")[1];
    repr := ReadRepr(filename);
    if repr = fail then
        return fail;
    fi;
    G := IndexSymmGroup(repr.m, repr.L, repr.ord);
    R := ExpandRepr(repr.M, G).L;
    write := WriteCDD(basename, repr.n, R, repr.t, repr.hrep);
    Print("Running scdd on ", write.cddinput, " ...\n");
    Exec(FileSCDD, " ", write.cddinput, " > /dev/null 2> /dev/null");
    read := ReadCDD(write.cddoutput);
    inc := ReadCDDInc(write.cddinc);
    M := CollapseRepr(read.R, G, inc);
    WriteRepr(basename, repr.m, repr.L, M, [], read.t, repr.ord, repr.diag, not repr.hrep);
    Print("Done.\n");
end;

CorrectSign := function (v, x, EXT, n)
    local j;
    j := Filtered([1..Length(EXT)], i -> not (i in x))[1];
    if EXT[j]{[2..n+1]} * v < 0 then
        return -v;
    else
        return v;
    fi;
end;

BuildGraph := function (G, inc, inc2)
    local i, j, P, O, p, o, n, L, A, Gr;
    Print("Entering BuildGraph...\n");
    n := Length(inc2);
    P := [];
    O := [];
    for i in [1..n] do
        for j in [i+1..n] do
            p := [i, j];
            if p in O then
                continue;
            fi;
            o := Orbit(G, p, OnSets);
            Add(P, rec(p := p, len := Length(o)));
            Append(O, o);
        od;
    od;
    L := Filtered(P, x -> Intersection(Union(List(Intersection(inc2[x.p[1]], inc2[x.p[2]]), y -> inc[y]), [[1..n]])) = x.p);
    A := Union(List(L, x -> Orbit(G, x.p, OnSets)));
    Gr := Graph(Group(()), [1..n], function (x, y) return x; end, function (x, y) return Set([x, y]) in A; end, true);
    return rec(adj := L, graph := Gr, adjmat := CollapsedAdjacencyMat(G, Gr), diam := Diameter(Gr));
end;

WriteGraph := function(basename, gd, m, L, t, g, hrep)
    local filename, f, l;
    
    if hrep then
        filename := Concatenation(basename, ".rdg");
    else
        filename := Concatenation(basename, ".skg");
    fi;
    Print("Writing file ", filename, " ...\n");
    f := OutputTextFile(filename, false);
    
    if hrep then
        WriteLine(f, "rdg_file: Ridge graph");
    else
        WriteLine(f, "skg_file: Skeleton graph");
    fi;
    WriteLine(f, Concatenation("name ", filename));
    
    WriteLine(f, Concatenation("vertices ", String(gd.graph.order)));
    WriteLine(f, Concatenation("diameter ", String(gd.diam)));
    
    WriteLine(f, "adjacencies");
    WriteLine(f, Concatenation("\t", String(Length(gd.adj)), " ", String(m), " ", t));
    l := JoinStringsWithSeparator(List(L,
        x -> JoinStringsWithSeparator(List(x, String), ",")), "\t");
    WriteLine(f, Concatenation("\t", l, "\t", l, "\ts"));
    
    Sort(gd.adj, SortByLen);
    for l in List(gd.adj, x -> [List(x.p, g), x.len]) do
        WriteLine(f, Concatenation("\t", JoinStringsWithSeparator(l[1][1], "\t"),
            "\t", JoinStringsWithSeparator(l[1][2], "\t"), "\t", String(l[2])));
    od;
    
    WriteLine(f, "end");
    
    WriteLine(f, "collapsed adjacency matrix");
    WriteLine(f, Concatenation("\t", String(Length(gd.adjmat))));
    for l in gd.adjmat do
        WriteLine(f, Concatenation("\t", JoinStringsWithSeparator(l, "\t")));
    od;
    WriteLine(f, "end");
    CloseStream(f);
end;

ADMConvert := function (filename)
    local basename, repr, G, exp, H, H2, f, g, EXT, nEXT, l, inc, inc2, conv, orb, M, N, res, res2, rd, sk, rh, sh;
    basename := SplitString(filename, ".")[1];
    repr := ReadRepr(filename);
    if repr = fail then
        return fail;
    fi;
    G := IndexSymmGroup(repr.m, repr.L, repr.ord);
    exp := ExpandRepr(repr.M, G);
    H := Action(G, exp.L, OnIndices);
    if repr.t = "integer" then
        f := Int;
        g := NullspaceIntMat;
    else
        f := Rat;
        g := NullspaceMat;
    fi;
    EXT := List(exp.L, x -> Concatenation([0], List(x, y -> f(y))));
    nEXT := NullspaceIntMat(TransposedMat(EXT{[1..Length(EXT)]}{[2..repr.n+1]}));
    inc := DualDescriptionStandard(EXT, H);
    conv := List(inc, x -> CorrectSign(g(TransposedMat(Union(EXT{x}{[2..repr.n+1]}, nEXT)))[1], x, EXT, repr.n));
    orb := List(inc, x -> Orbit(H, x, OnSets));
    res := Union(orb);
    res2 := List([1..Length(exp.L)], x -> Filtered([1..Length(res)], y -> x in res[y]));
    inc2 := res2{exp.pos};
    M := List([1..Length(inc)], x -> rec(repr := conv[x], inc := Length(inc[x]), len := Length(orb[x])));
    N := List([1..Length(exp.len)], x -> rec(repr := repr.M[x], inc := Length(inc2[x]), len := exp.len[x]));
    WriteRepr(basename, repr.m, repr.L, M, nEXT, repr.t, repr.ord, repr.diag, not repr.hrep);
    WriteRepr(basename, repr.m, repr.L, N, repr.N, repr.t, repr.ord, repr.diag, repr.hrep);
    Print("Done.\n");
    if repr.hrep then
        rd := res;
        sk := res2;
        rh := x -> Flat(EXT{[x]}{[2..repr.n+1]});
        sh := x -> CorrectSign(g(TransposedMat(Union(EXT{rd[x]}{[2..repr.n+1]}, nEXT)))[1], rd[x], EXT, repr.n);
    else
        rd := res2;
        sk := res;
        sh := x -> Flat(EXT{[x]}{[2..repr.n+1]});
        rh := x -> CorrectSign(g(TransposedMat(Union(EXT{sk[x]}{[2..repr.n+1]}, nEXT)))[1], sk[x], EXT, repr.n);
    fi;
    return rec(basename := basename, H := H, repr := repr, rd := rd, sk := sk,
        rh := rh, sh := sh, nullspace := nEXT);
end;

RidgeGraph := function (data)
    local H, gd;
    Print("Building ridge graph...\n");
    if data.repr.hrep then
        H := data.H;
    else
        H := Action(data.H, data.sk, OnSets);
    fi;
    gd := BuildGraph(H, data.rd, data.sk);
    WriteGraph(data.basename, gd, data.repr.m, data.repr.L, data.repr.t, data.rh, true);
end;

SkeletonGraph := function (data)
    local H, gd;
    Print("Building skeleton graph...\n");
    if data.repr.hrep then
        H := Action(data.H, data.rd, OnSets);
    else
        H := data.H;
    fi;
    gd := BuildGraph(H, data.sk, data.rd);
    WriteGraph(data.basename, gd, data.repr.m, data.repr.L, data.repr.t, data.sh, false);
end;


WPMET := function (m, basename)
    local L, R, i, j, k, l, z, r, len;
    L := MakeIndices(m, false, true);
    z := ListWithIdenticalEntries(Length(L), 0);
    i := Position(L, [1, 1]);
    R := [];
    r := ShallowCopy(z);
    r[i] := 1;
    Add(R, rec(repr := r, len := m));
    if m > 1 then
        r := ShallowCopy(z);
        r[i] := -1;
        j := Position(L, [1, 2]);
        if m = 2 then
            k := Position(L, [2, 2]);
            r[j] := 2;
            r[k] := -1;
            len := 1;
        else
            k := Position(L, [1, 3]);
            l := Position(L, [2, 3]);
            r[j] := 1;
            r[k] := 1;
            r[l] := -1;
            len := m*(m-1)*(m-2)/2;
        fi;
        Add(R, rec(repr := r, len := len));
    fi;
    
    WriteRepr(basename, m, L, R, [], "integer", false, true, true);
end;

PMET := function (m, basename)
    local L, R, i, j, k, l, z, r;
    L := MakeIndices(m, false, true);
    z := ListWithIdenticalEntries(Length(L), 0);
    i := Position(L, [1, 1]);
    R := [];
    r := ShallowCopy(z);
    r[i] := 1;
    Add(R, rec(repr := r, len := m));
    if m > 1 then
        r := ShallowCopy(z);
        j := Position(L, [1, 2]);
        r[i] := -1;
        r[j] := 1;
        Add(R, rec(repr := r));
        if m > 2 then
            r := ShallowCopy(r);
            k := Position(L, [1, 3]);
            l := Position(L, [2, 3]);
            r[k] := 1;
            r[l] := -1;
            Add(R, rec(repr := r, len := m*(m-1)*(m-2)/2));
        fi;
    fi;
    
    WriteRepr(basename, m, L, R, [], "integer", false, true, true);
end;

DWMET := function (m, basename)
    local L, R, i, j, k, l, z, r;
    L := MakeIndices(m, false, true);
    z := ListWithIdenticalEntries(Length(L), 0);
    i := Position(L, [1, 1]);
    R := [];
    r := ShallowCopy(z);
    r[i] := 1;
    Add(R, rec(repr := r, len := m));
    if m > 1 then
        r := ShallowCopy(z);
        j := Position(L, [1, 2]);
        k := Position(L, [2, 2]);
        r[i] := -1;
        r[j] := 1;
        r[k] := 1;
        Add(R, rec(repr := r, len := m*(m-1)));
        if m > 2 then
            r := ShallowCopy(z);
            k := Position(L, [1, 3]);
            l := Position(L, [2, 3]);
            r[j] := 1;
            r[k] := 1;
            r[l] := -1;
            Add(R, rec(repr := r, len := m*(m-1)*(m-2)/2));
        fi;
    fi;
    
    WriteRepr(basename, m, L, R, [], "integer", false, true, true);
end;

MET := function (m, basename)
    local L, R, i, j, k, r;
    L := MakeIndices(m, false, false);
    R := [];
    if m > 1 then
        r := ListWithIdenticalEntries(Length(L), 0);
        i := Position(L, [1, 2]);
        r[i] := 1;
        if m > 2 then
            j := Position(L, [1, 3]);
            k := Position(L, [2, 3]);
            r[j] := 1;
            r[k] := -1;
        fi;
        Add(R, rec(repr := r, len := m*(m-1)*(m-2)/2));
    fi;
    
    WriteRepr(basename, m, L, R, [], "integer", false, false, true);
end;

QMET := function (m, basename)
    local L, R, i, j, k, r;
    L := MakeIndices(m, true, false);
    R := [];
    if m > 1 then
        r := ListWithIdenticalEntries(Length(L), 0);
        i := Position(L, [1, 2]);
        r[i] := 1;
        Add(R, rec(repr := r, len := m*(m-1)));
        if m > 2 then
            r := ShallowCopy(r);
            j := Position(L, [1, 3]);
            k := Position(L, [3, 2]);
            r[j] := 1;
            r[k] := -1;
            Add(R, rec(repr := r, len := m*(m-1)*(m-2)));
        fi;
    fi;
    
    WriteRepr(basename, m, L, R, [], "integer", true, false, true);
end;

ZeroOneWMET := function (m, basename)
    local L, R, i, j, k, z, r, len;
    L := MakeIndices(m, false, true);
    R := [];
    z := ListWithIdenticalEntries(Length(L), 0);
    r := ShallowCopy(z);
    i := Position(L, [1, 1]);
    r[i] := 1;
    Add(R, rec(repr := r, len := m));
    for i in [1..Int(m/2)] do
        r := ShallowCopy(z);
        k := 0;
        for j in L do
            k := k+1;
            if j[1] <= i and j[2] > i then
                r[k] := 1;
            fi;
        od;
        len := Binomial(m, i);
        if i = m/2 then
            len := len/2;
        fi;
        Add(R, rec(repr := r, len := len));
    od;
    
    WriteRepr(basename, m, L, R, [], "integer", false, true, false);
end;

NMET := function (m, basename)
    local L, R, i, j, k, l, r, len;
    L := MakeIndices(m, false, true);
    r := ListWithIdenticalEntries(Length(L), 0);
    R := [];
    if m > 1 then
        i := Position(L, [1, 1]);
        j := Position(L, [1, 2]);
        r[i] := -1;
        if m = 2 then
            k := Position(L, [2, 2]);
            r[j] := 2;
            r[k] := -1;
            len := 1;
        else
            k := Position(L, [1, 3]);
            l := Position(L, [2, 3]);
            r[j] := 1;
            r[k] := 1;
            r[l] := -1;
            len := m*(m-1)*(m-2)/2;
        fi;
        Add(R, rec(repr := r, len := len));
    fi;
    
    WriteRepr(basename, m, L, R, [], "integer", false, true, true);
end;

NQMET := function (m, basename)
    local L, R, i, j, k, l, r, len;
    L := MakeIndices(m, true, true);
    r := ListWithIdenticalEntries(Length(L), 0);
    R := [];
    if m > 1 then
        i := Position(L, [1, 1]);
        j := Position(L, [2, 1]);
        if m = 2 then
            k := Position(L, [1, 2]);
            l := Position(L, [2, 2]);
            len := 1;
        else
            k := Position(L, [1, 3]);
            l := Position(L, [2, 3]);
            len := m*(m-1)*(m-2);
        fi;
        r[i] := -1;
        r[j] := 1;
        r[k] := 1;
        r[l] := -1;
        Add(R, rec(repr := r, len := len));
    fi;
    
    WriteRepr(basename, m, L, R, [], "integer", true, true, true);
end;

CUT := function (m, basename)
    local L, R, i, j, k, r, len;
    L := MakeIndices(m, false, false);
    R := [];
    for i in [1..Int(m/2)] do
        r := ListWithIdenticalEntries(Length(L), 0);
        k := 0;
        for j in L do
            k := k+1;
            if j[1] <= i and j[2] > i then
                r[k] := 1;
            fi;
        od;
        len := Binomial(m, i);
        if i = m/2 then
            len := len/2;
        fi;
        Add(R, rec(repr := r, len := len));
    od;
    
    WriteRepr(basename, m, L, R, [], "integer", false, false, false);
end;

OCUT := function (m, basename)
    local L, R, i, j, k, r;
    L := MakeIndices(m, true, false);
    R := [];
    for i in [1..m-1] do
        r := ListWithIdenticalEntries(Length(L), 0);
        k := 0;
        for j in L do
            k := k+1;
            if j[1] <= i and j[2] > i then
                r[k] := 1;
            fi;
        od;
        Add(R, rec(repr := r, len := Binomial(m, i)));
    od;
    
    WriteRepr(basename, m, L, R, [], "integer", true, false, false);
end;

ZeroOne := function (filename)
    local basename, repr, f, M;
    basename := SplitString(filename, ".")[1];
    repr := ReadRepr(filename);
    if repr = fail then
        return fail;
    fi;
    if repr.t = "integer" then
        f := Int;
    else
        f := Rat;
    fi;
    M := List(Filtered(repr.M,
        x -> Length(Filtered(x, y -> not (f(y) in [0, 1]))) = 0),
        z -> rec(repr := z));
    WriteRepr(Concatenation("01-", basename), repr.m, repr.L, M, [], repr.t, repr.ord, repr.diag, repr.hrep);
    Print("Done.\n");
end;