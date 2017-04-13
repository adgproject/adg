%%
%% This is file `feynmp.mp',
%% generated with the docstrip utility.
%%
%% The original source files were:
%%
%% feynmf.dtx  (with options: `base,mp')
%% 
%% Copyright (C) 1989, 1990, 1992-1995 by Thorsten.Ohl@Physik.TH-Darmstadt.de
%% 
%% This file is NOT the source for feynmf, because almost all comments
%% have been stripped from it.  It is NOT the preferred form of feynmf
%% for making modifications to it.
%% 
%% Therefore you can NOT redistribute and/or modify THIS file.  You can
%% however redistribute the complete source (feynmf.dtx and feynmf.ins)
%% and/or modify it under the terms of the GNU General Public License as
%% published by  the Free Software Foundation; either version 2, or (at
%% your option) any later version.
%% 
%% As a special exception, you can redistribute parts of this file for
%% the electronic distribution of scientific papers, provided that you
%% include a short note pointing to the complete source.
%% 
%% Feynmf is distributed in the hope that it will be useful, but
%% WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%% 
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% \CheckSum{924}
%% \CharacterTable
%%  {Upper-case    \A\B\C\D\E\F\G\H\I\J\K\L\M\N\O\P\Q\R\S\T\U\V\W\X\Y\Z
%%   Lower-case    \a\b\c\d\e\f\g\h\i\j\k\l\m\n\o\p\q\r\s\t\u\v\w\x\y\z
%%   Digits        \0\1\2\3\4\5\6\7\8\9
%%   Exclamation   \!     Double quote  \"     Hash (number) \#
%%   Dollar        \$     Percent       \%     Ampersand     \&
%%   Acute accent  \'     Left paren    \(     Right paren   \)
%%   Asterisk      \*     Plus          \+     Comma         \,
%%   Minus         \-     Point         \.     Solidus       \/
%%   Colon         \:     Semicolon     \;     Less than     \<
%%   Equals        \=     Greater than  \>     Question mark \?
%%   Commercial at \@     Left bracket  \[     Backslash     \\
%%   Right bracket \]     Circumflex    \^     Underscore    \_
%%   Grave accent  \`     Left brace    \{     Vertical bar  \|
%%   Right brace   \}     Tilde         \~}
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def subgraph (expr x, y, wd, ht) =
  begingroup
    save c, ne, nw, sw, se;
    pair c, ne, nw, sw, se;
    sw = (x,y);
    se = sw + (wd,0);
    nw = sw + (0,ht);
    ne = sw + (wd,ht);
    c = .5[sw,ne];
enddef;
def endsubgraph =
  endgroup
enddef;
if known cmbase:
  errhelp
    "feynmf will only work with plain Metafont, as described in the book.";
  errmessage "feynmf: CMBASE detected.  Please use the PLAIN base.";
  forever:
    errmessage "No use in trying!  You'd better eXit now ...";
    errorstopmode;
  endfor
fi
vardef parse_RCS (suffix RCS) (expr s) =
  save n, c;
  numeric n, RCS[];
  string c;
  RCS[0] := 0;
  for n = 1 upto length (s):
    c := substring (n-1,n) of s;
    exitif ((RCS[0] > 0) and (c = " "));
    if ((c = "0") or (c = "1") or (c = "2")
        or (c = "3") or (c = "4") or (c = "5")
        or (c = "6") or (c = "7") or (c = "8")
        or (c = "9")):
      if RCS[0] = 0:
        RCS[0] := 1;
        RCS[RCS[0]] := 0;
      fi
      RCS[RCS[0]] := 10 * RCS[RCS[0]] + scantokens (c);
    elseif c = ".":
      RCS[0] := RCS[0] + 1;
      RCS[RCS[0]] := 0;
    else:
    fi
  endfor
enddef;
vardef require_RCS_revision expr s =
  save n, TeX_rev, mf_rev;
  numeric n;
  parse_RCS (TeX_rev, s);
  parse_RCS (mf_rev, "$Revision: 1.30 $");
  for n = 1 upto min (2, TeX_rev[0], mf_rev[0]):
    if TeX_rev[n] > mf_rev[n]:
      errhelp
        "Your version of `feynmf.sty' is higher that of your `feynmf.mf'.";
      errmessage "feynmf: Metafont macros out of date";
    elseif TeX_rev[n] < mf_rev[n]:
      errhelp
        "Your version of `feynmf.mf' is higher that of your `feynmf.sty'.";
      errmessage "feynmf: LaTeX style out of date";
    fi
    exitif (TeX_rev[n] <> mf_rev[n]);
  endfor
enddef;
vardef cullit = \ enddef;
color foreground;
foreground = black;
vardef beginchar (expr c, wd, ht, dp) =
  LaTeX_file := "";
  beginfig(c);
    w:=wd;
    h:=ht;
enddef;
string LaTeX_file;
vardef endchar =
  setbounds currentpicture to (0,0)--(w,0)--(w,h)--(0,h)--cycle;
  if LaTeX_file <> "":
    write EOF to LaTeX_file;
    LaTeX_file := "";
  fi
  endfig
enddef;
bp# := bp;
cc# := cc;
cm# := cm;
dd# := dd;
in# := in;
mm# := mm;
pc# := pc;
pt# := pt;
vardef define_blacker_pixels(text t) =
  forsuffixes $=t:
    $:=$.#;
  endfor
enddef;
picture unitpixel;
unitpixel = nullpicture;
addto unitpixel contour unitsquare;
def t_ = \ enddef;
boolean feynmfwizard;
feynmfwizard := false;
thin# := 1pt#; % dimension of the lines
thick# := 2thin#;
arrow_len# := 4mm#;
arrow_ang := 15;
curly_len# := 3mm#;
dash_len# := 3mm#; % 'photon' lines
dot_len# := 2mm#; % 'photon' lines
wiggly_len# := 4mm#; % 'photon' lines
wiggly_slope := 60;
zigzag_len# := 2mm#;
zigzag_width# := 2thick#;
decor_size# := 5mm#;
dot_size# := 2thick#;
define_blacker_pixels (thick, thin,
  dash_len, dot_len, wiggly_len, curly_len,
  zigzag_len, zigzag_width, arrow_len, decor_size, dot_size);
def shrink expr s =
  begingroup
  if shrinkables <> "":
    save tmp_;
    forsuffixes $ = scantokens shrinkables:
      if known $.#:
        tmp_ := $.#;
        save $;
        $.# := s * tmp_;
        define_blacker_pixels ($);
      else:
        tmp_ := $;
        save $;
        $ := s * tmp_;
      fi
    endfor
  fi
enddef;
def endshrink =
  endgroup
enddef;
string shrinkables;
shrinkables := "";
vardef addto_shrinkables (text l) =
  forsuffixes $ = l:
    shrinkables := shrinkables & "," & str $;
  endfor
enddef;
shrinkables := "thick,thin";
addto_shrinkables (dash_len, dot_len);
addto_shrinkables (wiggly_len, curly_len);
addto_shrinkables (zigzag_len, zigzag_width);
addto_shrinkables (arrow_len);
addto_shrinkables (decor_size, dot_size);
LaTeX_unitlength := mm;
vardef count (text list) =
  forsuffixes $ = list: + 1 endfor
enddef;
vardef getopt (suffix opt) (expr s) =
  save n, argp, escape, anchor, skip;
  numeric opt.first, opt.last, n, anchor;
  string opt[], opt[]arg;
  boolean opt[]tainted, argp, escape, skip;
  opt.first := 0;
  opt.last := 0;
  opt[opt.last] := "";
  argp := false;
  escape := false;
  anchor := 0;
  skip := true;
  for n = 1 upto length (s):
    if skip and (substring (n-1, n) of s = " "):
      anchor := anchor + 1;
    else:
      skip := false;
      if not escape and (substring (n-1, n) of s = ","):
        if substring (n, n+1) of s = ",":
          escape := true;
          opt[opt.last]tainted := true;
        else:
          if argp:
            opt[opt.last]arg := substring (anchor, n-1) of s;
          else:
            opt[opt.last] := substring (anchor, n-1) of s;
          fi
          anchor := n;
          argp := false;
          skip := true;
          opt.last := opt.last + 1;
        fi
      elseif not argp and (substring (n-1, n) of s = "="):
        opt[opt.last] := substring (anchor, n-1) of s;
        anchor := n;
        argp := true;
        skip := true;
      elseif argp or (substring (n-1, n) of s <> " "):
        escape := false;
      fi
    fi
  endfor
  if argp:
    opt[opt.last]arg := substring (anchor, length s) of s;
  else:
    opt[opt.last] := substring (anchor, length s) of s;
  fi
  for n = opt.first upto opt.last:
    if known opt[n]tainted:
      if opt[n]tainted:
        opt[n]arg := untaint_string opt[n]arg;
      fi
    fi
  endfor
enddef;
vardef untaint_string suffix s =
  save n, anchor;
  numeric n, anchor;
  anchor := 0;
  for n = 1 upto length (s) - 1:
    if substring (n-1,n+1) of s = ",,":
      substring (anchor, n-1) of s &
      hide (anchor := n)
    fi
  endfor
  substring (anchor, length s) of s
enddef;
vardef split_string (suffix comp) (expr s) =
  save n, anchor;
  numeric comp.first, comp.last, n, anchor;
  string comp[];
  comp.first := 0;
  comp.last := 0;
  comp[comp.last] := "";
  anchor := 0;
  for n = 1 upto length (s):
    if substring (n-1,n) of s = ".":
      comp[comp.last] := substring (anchor, n-1) of s;
      comp.last := comp.last + 1;
      anchor := n;
    fi
  endfor
  comp[comp.last] := substring (anchor, length s) of s;
enddef;
vardef match_prefix (expr prefix, s) =
  (prefix = substring (0, length prefix) of s)
enddef;
vardef match_option (expr s, option) =
  save sc, optionc, n, i;
  numeric sc.first, sc.last, optionc.first, optionc.last;
  string sc[], optionc[];
  numeric n, i;
  split_string (sc, s);
  split_string (optionc, option);
  n := sc.last - sc.first;
  if n <> (optionc.last - optionc.first):
    false
  else:
    true
    for i = 0 upto n:
      and match_prefix (sc[sc.first+i],
                        optionc[optionc.first+i])
    endfor
  fi
enddef;
def save_picture text t =
 save t; picture t; forsuffixes p=t: p:=nullpicture; endfor
enddef;
def begin_sketch =
 begingroup save_picture currentpicture;
 sketchlevel := sketchlevel+1;
enddef;
def end_sketch =
 sketchlevel := sketchlevel-1;
 sketchpad[sketchlevel] := currentpicture;
 endgroup
enddef;
picture sketchpad[];
sketchlevel := 1;
vardef use_sketch text t =
 addto currentpicture also (sketchpad[sketchlevel] t)
enddef;
vardef cdraw expr p =
  draw p withcolor foreground
enddef;
vardef cfill expr p =
  fill p withcolor foreground
enddef;
vardef cfilldraw expr p =
  filldraw p withcolor foreground
enddef;
vardef ccutdraw expr p =
  cutdraw p withcolor foreground
enddef;
vardef cdrawdot expr p =
  drawdot p withcolor foreground
enddef;
vardef isdigit expr s =
  save n;
  (s = "0")
  for n = 1 upto 9:
    or (s = decimal n)
  endfor
enddef;
vardef digits_index (expr s, start) =
  save n, m, from, to;
  for n = start upto (length s)-1:
    if isdigit (substring (n,n+1) of s):
      from := n;
      for m = n upto length s:
        if not isdigit (substring (m,m+1) of s):
          to := m;
        fi
        exitif known to;
      endfor
    fi
    exitif known from;
  endfor
  (from, if known to: to else: infinity fi)
enddef;
vardef digits_to_brackets suffix suf =
  save s, idx;
  string s;
  pair idx;
  s = str suf;
  idx = (0,0);
  forever:
    idx := digits_index (s, xpart idx);
    exitif unknown xpart idx;
    s := substring (0,xpart idx) of s
      & "[]" & substring (ypart idx,infinity) of s;
  endfor
  s
enddef;
tile_grain := 1in/300;
vardef def_tile (suffix t) (expr wd, ht) =
  if not picture tlist.t:
    picture tlist.scantokens (digits_to_brackets t);
  fi
  tlist.t := nullpicture;
  tlist.t.dx := max (floor wd, 1);
  tlist.t.dy := max (floor ht, 1);
enddef;
vardef use_tile (suffix t) (expr x, y, wd, ht) =
  fill unitsquare xscaled wd yscaled ht shifted (x,y)
    withcolor background;
  if str t = "shaded":
    shade_rectangle (4thin, x, y, wd, ht);
  elseif str t = "hatched":
    shade_rectangle (5thin, x, y, wd, ht);
    shade_rectangle (-5thin, x, y, wd, ht);
  else:
    if (picture tlist.t):
      for nx = 0 upto wd/tlist.t.dx:
        for ny = 0 upto ht/tlist.t.dy:
          addto currentpicture
            also (tlist.t shifted
                  ((x,y) + (nx*tlist.t.dx, ny*tlist.t.dy)) t_);
        endfor
      endfor
    else:
      errhelp "feynmf: your tiling has not been defined, "
            & "check spelling and reprocess!";
      errmessage "feynmf: tiling `" & str t & "' not known, "
               & "replaced by `shaded'";
      use_tile (shaded, x, y, wd, ht);
    fi
  fi
enddef;
vardef shade_rectangle (expr dd, x, y, wd, ht) =
  save d, u, dx, dy, currentpen;
  pen currentpen;
  pickup pencircle scaled thin;
  d := max (floor (abs dd), 1);
  dx := max (wd, ht);
  dy := max (wd, ht);
  for u = 0 step d/dx until 1:
    if dd > 0:
      cdraw (x-d,y+u*dy-d)--(x+(1-u)*dx+d,y+dy+d);
      cdraw (x+u*dx-d,y-d)--(x+dx+d,y+(1-u)*dy+d);
    else:
      cdraw (x-d,y+u*dy+d)--(x+u*dx+d,y-d);
      cdraw (x+(1-u)*dx-d,y+dy+d)--(x+dx+d,y+(1-u)*dy-d);
    fi
  endfor
enddef;
def addto_tile (suffix t) =
  addto tlist.t
enddef;
vardef tile_from_string (suffix t) (expr str) =
  tile_grain := max (floor tile_grain, 1);
  save grain, mx, x, y, n, c, pic;
  string c;
  picture pic;
  pic := nullpicture;
  grain := tile_grain;
  mx := 0;
  x := 0;
  y := 0;
  for n := 1 upto length str:
    c := substring (n-1,n) of str;
    if c = "/":
      mx := max (mx, x);
      y := y+1;
      x := 0;
    elseif c = "*":
      addto pic also (unitpixel shifted (x,-y) t_);
      x := x+1;
    elseif c = ".":
      x := x+1;
    fi
  endfor
  def_tile (t, grain*mx, grain*(y+1));
  addto_tile (t) also (pic shifted (0,y) scaled grain t_);
  pic := nullpicture;
enddef;
tile_from_string (gray10,
  "  ...  /"&
  "  .*.  /"&
  "  ...  ");
tile_from_string (gray25,
  "  ..  /"&
  "  *.  ");
tile_from_string (gray50,
  "  .*  /"&
  "  *.  ");
tile_from_string (gray75,
  "  **  /"&
  "  .*  ");
tile_from_string (gray90,
  "  ***  /"&
  "  *.*  /"&
  "  ***  ");
vardef make_halftone (suffix t) (expr g, wd, ht) =
  def_tile (t, wd, ht);
  addto tlist.t contour unitsquare
    xscaled wd yscaled ht withcolor ((1-g)*foreground + g*background)
enddef;
vardef tile (suffix t) (expr p) =
  save u, x, y, max_x, min_x, max_y, min_y, xx, yy;
  -max_x = -max_y = min_x = min_y = infinity;
  for u = 0 step 0.1 until length p:
    x := xpart (point u of p);
    y := ypart (point u of p);
    max_x := max(max_x, x);
    max_y := max(max_y, y);
    min_x := min(min_x, x);
    min_y := min(min_y, y);
  endfor
  begin_sketch
    use_tile (t, min_x, min_y, max_x-min_x, max_y-min_y);
    clip currentpicture to p;
  end_sketch;
  use_sketch;
enddef;
vardef drawtile (suffix t) (expr p) =
  tile (t, p);
  cdraw p
enddef;
vardef use_halftone (expr g, x, y, wd, ht) =
  fill unitsquare xscaled wd yscaled ht
    shifted (x,y) withcolor ((1-g)*foreground + g*background)
enddef;
vardef halftone (expr g, p) =
  fill p withcolor ((1-g)*foreground + g*background)
enddef;
vardef drawhalftone (expr g, p) =
  fill p withcolor ((1-g)*foreground + g*background);
  cdraw p
enddef;
vardef shade expr p =
  tile (shaded, p)
enddef;
vardef hatch expr p =
  tile (hatched, p)
enddef;
vardef emptydraw expr p =
  cullit;
  unfill p;
  cullit;
  cdraw p;
enddef;
vardef shadedraw expr p =
  cullit;
  unfill p;
  cullit;
  shade p;
  cdraw p;
enddef;
vardef hatchdraw expr p =
  cullit;
  unfill p;
  cullit;
  hatch p;
  cdraw p;
enddef;
vardef marrow (expr p, frac) =
  save a, t, z;
  pair z;
  a = angle direction frac*length(p) of p;
  z = point frac*length(p) of p;
  (t1,whatever) = p intersectiontimes
    (halfcircle scaled 2/3arrow_len rotated (a+90) shifted z);
  (t2,whatever) = p intersectiontimes
    (halfcircle scaled 4/3arrow_len rotated (a-90) shifted z);
  arrow_head (p, t1, t2, arrow_ang)
enddef;
vardef tarrow (expr p, frac) =
  save a, t, z;
  pair z;
  t1 = frac*length p;
  a = angle direction t1 of p;
  z = point t1 of p;
  (t2,whatever) = p intersectiontimes
    (halfcircle scaled 2arrow_len rotated (a-90) shifted z);
  arrow_head (p, t1, t2, arrow_ang)
enddef;
vardef harrow (expr p, frac) =
  save a, t, z;
  pair z;
  t2 = frac*length p;
  a = angle direction t2 of p;
  z = point t2 of p;
  (t1,whatever) = p intersectiontimes
     (halfcircle scaled 2arrow_len rotated (a+90) shifted z);
  arrow_head (p, t1, t2, arrow_ang)
enddef;
vardef arrow_head (expr p, from, to, ang) =
  save tip, ap, t;
  pair tip;
  path ap;
  t1 := from;
  t2 := to;
  if t1 = -1: t1 := 0; fi
  if t2 = -1: t2 := infinity; fi
  tip = point t2 of p;
  ap = subpath (t1,t2) of p shifted -tip;
  (ap rotated ang
    forced_join reverse ap rotated -ang
    -- cycle) shifted tip
enddef;
vardef arrow expr p =
  marrow (p, .5)
enddef;
tertiarydef p forced_join q =
  subpath (0, length p - 1) of p
  & point (length p - 1) of p
    .. controls postcontrol (length p - 1) of p
                and precontrol infinity of p
  .. .5[point infinity of p, point 0 of q]
    .. controls postcontrol 0 of q and precontrol 1 of q
    .. point 1 of q
  & subpath (1, infinity) of q
enddef;
vardef cut_decors (suffix from) (expr p) (suffix to) =
 subpath (if known from.decor.shape:
            xpart (p intersectiontimes
                     (from.decor.shape scaled from.decor.size
                                       shifted from.loc))
          else:
            0
          fi,
          if known to.decor.shape:
            length p
              - xpart (reverse p intersectiontimes
                        (to.decor.shape scaled to.decor.size
                                        shifted to.loc))
          else:
             infinity
          fi) of p
enddef;
vardef make_blob (expr z_arg, diameter) =
 save p,currentpen; path p; pen currentpen;
 pickup pencircle scaled thick;
 p = fullcircle scaled diameter shifted z_arg;
 shadedraw p;
enddef;
vardef draw_blob (expr z_arg, diameter) =
 if sketched_blob_diameter <> diameter: % drawn lately?
  begin_sketch make_blob (origin, diameter); end_sketch; % redo hard work!
  sketched_blob_diameter:= diameter;  % record it
 fi
 use_sketch shifted z_arg; % the easy way ...
enddef;
def force_new_blob = sketched_blob_diameter := -1; enddef;
force_new_blob;                                 % initialize it.
vardef pixlen (expr p, n) =
  for k=1 upto length(p): + segment_pixlen (subpath (k-1,k) of p, n) endfor
enddef;
vardef segment_pixlen (expr p, n) =
  for k=1 upto n: + abs (point k/n of p - point (k-1)/n of p) endfor
enddef;
vardef wiggly expr p_arg =
 save wpp;
 numeric wpp;
 wpp = ceiling (pixlen (p_arg, 10) / wiggly_len) / length p_arg;
 for k=0 upto wpp*length(p_arg) - 1:
  point k/wpp of p_arg
       {direction k/wpp of p_arg rotated wiggly_slope} ..
  point (k+.5)/wpp of p_arg
       {direction (k+.5)/wpp of p_arg rotated - wiggly_slope} ..
 endfor
 if cycle p_arg: cycle else: point infinity of p_arg fi
enddef;
vardef curly expr p =
 save cpp;
 numeric cpp;
 cpp := ceiling (pixlen (p, 10) / curly_len) / length p;
 if cycle p:
   for k=0 upto cpp*length(p) - 1:
     point (k+.33)/cpp of p
           {direction (k+.33)/cpp of p rotated 90} ..
     point (k-.33)/cpp of p
           {direction (k-.33)/cpp of p rotated -90} ..
   endfor
   cycle
 else:
   point 0 of p
         {direction 0 of p rotated -90} ..
   for k=1 upto cpp*length(p) - 1:
     point (k+.33)/cpp of p
           {direction (k+.33)/cpp of p rotated 90} ..
     point (k-.33)/cpp of p
           {direction (k-.33)/cpp of p rotated -90} ..
   endfor
   point infinity of p
         {direction infinity of p rotated 90}
 fi
enddef;
vardef zigzag expr p =
 save zpp;
 numeric zpp;
 zpp = ceiling (pixlen (p, 10) / zigzag_len) / length p;
 if not cycle p:
   point 0 of p --
 fi
 for k = 0 upto zpp*length(p) - 1:
   point (k+1/3)/zpp of p shifted
     (zigzag_width
      * dir angle (direction (k+1/3)/zpp of p rotated 90)) --
   point (k+2/3)/zpp of p shifted
     (zigzag_width
      * dir angle (direction (k+2/3)/zpp of p rotated -90)) --
 endfor
 if cycle p:
   cycle
 else:
   point infinity of p
 fi
enddef;
save vsty_hash;
def style_def suffix s =
  vsty_hash.s := 1;
  expandafter quote vardef scantokens ("draw_" & str s)
enddef;
vardef vsty_exists suffix s =
  known vsty_hash.s
enddef;
vardef valid_style expr s =
  expandafter vsty_exists scantokens (s)
enddef;
style_def phantom expr p =
  \
enddef;
style_def phantom_arrow expr p =
  cfill (arrow p);
enddef;
style_def plain expr p =
  cdraw p;
enddef;
style_def plain_arrow expr p =
  cdraw p;
  cfill (arrow p);
enddef;
style_def dbl_plain expr p =
  draw_double p;
enddef;
style_def dbl_plain_arrow expr p =
  draw_double_arrow p;
enddef;
style_def wiggly expr p =
  cdraw (wiggly p);
enddef;
style_def dbl_wiggly expr p =
  draw_double (wiggly p);
enddef;
style_def curly expr p =
  cdraw (curly p);
enddef;
style_def dbl_curly expr p =
  draw_double (curly p);
enddef;
style_def zigzag expr p =
  cdraw (zigzag p);
enddef;
style_def dbl_zigzag expr p =
  draw_double (zigzag p);
enddef;
style_def dashes expr p =
 save dpp;
 numeric dpp;
 dpp = ceiling (pixlen (p, 10) / dash_len) / length p;
 for k=0 upto dpp*length(p) - 1:
  cdraw point k/dpp of p ..
   point (k+.5)/dpp of p;
 endfor
enddef;
style_def dbl_dashes expr p =
 save dpp;
 numeric dpp;
 dpp = ceiling (pixlen (p, 10) / dash_len) / length p;
 for k=0 upto dpp*length(p) - 1:
  draw_double point k/dpp of p ..
   point (k+.5)/dpp of p;
 endfor
enddef;
style_def dbl_dashes_arrow expr p =
  draw_dbl_dashes p;
  shrink (1.5);
    cfill (arrow p);
  endshrink;
enddef;
style_def dashes_arrow expr p =
  draw_dashes p;
  cfill (arrow p);
enddef;
style_def dots expr p =
 save dpp;
 numeric dpp;
 dpp = ceiling (pixlen (p, 10) / dot_len) / length p;
 for k=0 upto dpp*length(p):
   cdrawdot point k/dpp of p;
 endfor
enddef;
style_def dbl_dots expr p =
  save dpp;
  numeric dpp;
  dpp = ceiling (pixlen (p, 10) / dot_len) / length p;
  begingroup
    save oldpen;
    pen oldpen;
    oldpen := currentpen;
    pickup oldpen scaled 3; % draw a thick linn
    for k=0 upto dpp*length(p):
      cdrawdot point k/dpp of p;
    endfor
    pickup oldpen;
    cullit;
    for k=0 upto dpp*length(p):
      undrawdot point k/dpp of p;
    endfor
    cullit; % and remove the stuffing
  endgroup;
enddef;
style_def dbl_dots_arrow expr p =
  draw_dbl_dots p;
  shrink (1.5);
    cfill (arrow p);
  endshrink;
enddef;
style_def dots_arrow expr p =
  draw_dots p;
  cfill (arrow p);
enddef;
style_def double expr p =
  save oldpen;
  pen oldpen;
  oldpen := currentpen;
  pickup oldpen scaled 3;
  ccutdraw p;
  pickup oldpen;
  cullit; undraw p; cullit;
enddef;
style_def double_arrow expr p =
  draw_double p;
  shrink (1.5);
    cfill (arrow p);
  endshrink;
enddef;
style_def vanilla expr p = draw_plain p enddef;
style_def fermion expr p = draw_plain_arrow p enddef;
style_def quark expr p = draw_plain_arrow p enddef;
style_def electron expr p = draw_plain_arrow p enddef;
style_def photon expr p = draw_wiggly p enddef;
style_def boson expr p = draw_wiggly p enddef;
style_def gluon expr p = draw_curly p enddef;
style_def heavy expr p = draw_dbl_plain_arrow p enddef;
style_def ghost expr p = draw_dots_arrow p enddef;
style_def scalar expr p = draw_dashes_arrow p enddef;
% Homemade routines for graphs
style_def prop_pm expr p =
    draw_plain p;
    shrink(.7);
        cfill (marrow (p, .25));
        cfill (marrow (p, .75))
    endshrink;
enddef;
%
style_def prop_pp expr p =
    cdraw p;
    shrink(.7);
        cfill (marrow (p, .25));
        cfill (marrow (reverse p, .25))
    endshrink;
enddef;
%
style_def prop_mm expr p =
    draw_plain p;
    shrink(.7);
        cfill (marrow (p, .75));
        cfill (marrow (reverse p, .75))
    endshrink;
enddef;
%
style_def prop_mp expr p =
    draw_plain p;
    shrink(.7);
        cfill (marrow (reverse p, .25));
        cfill (marrow (reverse p, .75))
    endshrink;
enddef;
%
style_def half_prop expr p =
    draw_plain p;
    shrink(.7);
        cfill (marrow (p, .5))
    endshrink;
enddef;
%
vardef fermion expr path_arg =
  cfill (arrow (path_arg));
  path_arg
enddef;
vardef photon expr path_arg =
  wiggly path_arg
enddef;
vardef gluon expr path_arg =
  curly path_arg
enddef;
tracingstats:=1;
boolean vtracing;
vtracing := false; % true
def vinit =
  save vhash;
  numeric vlist.first, vlist.last;
  vlist.first := 1;
  vlist.last := 0;
  pair vlist[]loc;
  numeric vlist[]decor.size, vlist[]decor.ang,
    vlist[]arc.first, vlist[]arc.last,
    vlist[]arc[], vlist[]arc[]lsr,
    vlist[]arc[]tns, vlist[]arc[]lbl.dist,
    vlist[]arc[]tag, vlist[]arc[]wd, vlist[]arc[]rub,
    vlist[]constr.first, vlist[]constr.last,
    vlist[]constr[], lambdax[][], lambday[][];
  string vlist[]name, vlist[]lbl, vlist[]decor.sty,
    vlist[]arc[]sty, vlist[]arc[]lbl, vlist[]arc[]lbl.side;
  numeric vlist[]lbl.ang;
  path vlist[]decor.shape;
  color vlist[]fore, vlist[]back,
    vlist[]arc[]fore, vlist[]arc[]back;
  numeric plist.first, plist.last, plist[]cnt, plist[]vtx[],
   plist[]pull, plist[]lbl.ang, plist[]lbl.dist;
  string plist[]lbl, plist[]sty, plist[]cona, plist[]conb;
  plist.first := 1;
  plist.last := 0;
  numeric vlist[]poly.first, vlist[]poly.last,
    vlist[]poly[], vlist[]poly[]idx;
  pair lambdap[][];
 color plist[]fore, plist[]back;
enddef;
def vertices =
  vlist.first upto vlist.last
enddef;
def varcs (text i) =
  vlist[i]arc.first upto vlist[i]arc.last
enddef;
def vconstr (text i) =
  vlist[i]constr.first upto vlist[i]constr.last
enddef;
def polygons =
  plist.first upto plist.last
enddef;
def vpoly (text i) =
  vlist[i]poly.first upto vlist[i]poly.last
enddef;
vardef venter suffix v =
  if not vexists v:
    vlist.last := vlist.last + 1;
    vhash.v := vlist.last;
    vlist[vhash.v]name := str v;
    vlist[vhash.v]loc := (whatever,whatever);
    vlist[vhash.v]arc.first := 1;
    vlist[vhash.v]arc.last := 0;
    vlist[vhash.v]constr.first := 1;
    vlist[vhash.v]constr.last := 0;
    vlist[vhash.v]lbl := "";
    vlist[vhash.v]lbl.ang := whatever;
    vlist[vhash.v]lbl.dist := 3thick;
    vlist[vhash.v]fore := (whatever, whatever, whatever);
    vlist[vhash.v]back := (whatever, whatever, whatever);
    vlist[vhash.v]poly.first := 1;
    vlist[vhash.v]poly.last := 0;
  fi
enddef;
vardef vexists suffix v =
  if known vhash.v: true else: false fi
enddef;
vardef vlookup suffix v =
  if vexists v: vhash.v else: 0 fi
enddef;
vardef vloc suffix v =
  vlist[vlookup v]loc
enddef;
vardef vconnect (expr linesty) (text vl) =
  save from, nfrom, nto, nopt, sty;
  numeric from, nfrom, nto, nopt;
  string sty;
  getopt (opt, linesty);
  sty := opt[opt.first];
  if known opt[opt.first]arg:
    message "feynmf: line styles don't take arguments.  "
             & "Argument `" & opt[opt.first]arg & "' ignored.";
  fi
  opt.first := opt.first + 1;
  forsuffixes to = vl:
    venter to;
    nto := vlookup to;
    if known nfrom:
      vlist[nfrom]arc.last := vlist[nfrom]arc.last + 1;
      vlist[nfrom]arc[vlist[nfrom]arc.last] := nto;
      vlist[nfrom]arc[vlist[nfrom]arc.last]tns := 1;
      if nfrom <> nto:
        vlist[nto]arc.last := vlist[nto]arc.last + 1;
        vlist[nto]arc[vlist[nto]arc.last] := nfrom;
        vlist[nto]arc[vlist[nto]arc.last]tns := 1;
      fi
      vlist[nfrom]arc[vlist[nfrom]arc.last]lbl := "";
      vlist[nfrom]arc[vlist[nfrom]arc.last]lbl.side := "";
      vlist[nfrom]arc[vlist[nfrom]arc.last]lbl.dist := 3thick;
      for nopt = opt.first upto opt.last:
        if match_option (opt[nopt], "tension"):
          get_argument (opt[nopt], scantokens (opt[nopt]arg),
                        vlist[nfrom]arc[vlist[nfrom]arc.last]tns);
          get_argument (opt[nopt], scantokens (opt[nopt]arg),
                        vlist[nto]arc[vlist[nto]arc.last]tns);
        elseif match_option (opt[nopt], "left"):
          if known opt[nopt]arg:
            vlist[nfrom]arc[vlist[nfrom]arc.last]lsr
              := - scantokens (opt[nopt]arg);
          else:
            vlist[nfrom]arc[vlist[nfrom]arc.last]lsr := -1;
          fi
        elseif match_option (opt[nopt], "straight"):
          vlist[nfrom]arc[vlist[nfrom]arc.last]lsr := 0;
          ignore_argument (opt[nopt], opt[nopt]arg);
        elseif match_option (opt[nopt], "right"):
          if known opt[nopt]arg:
            vlist[nfrom]arc[vlist[nfrom]arc.last]lsr
              := scantokens (opt[nopt]arg);
          else:
            vlist[nfrom]arc[vlist[nfrom]arc.last]lsr := 1;
          fi
        elseif match_option (opt[nopt], "label"):
          get_argument (opt[nopt], opt[nopt]arg,
                        vlist[nfrom]arc[vlist[nfrom]arc.last]lbl);
        elseif match_option (opt[nopt], "label.side"):
          get_argument (opt[nopt], opt[nopt]arg,
                        vlist[nfrom]arc[vlist[nfrom]arc.last]lbl.side);
        elseif match_option (opt[nopt], "label.dist"):
          get_argument (opt[nopt], scantokens (opt[nopt]arg),
                        vlist[nfrom]arc[vlist[nfrom]arc.last]lbl.dist);
        elseif match_option (opt[nopt], "tag"):
          if known opt[nopt]arg:
            vlist[nfrom]arc[vlist[nfrom]arc.last]tag
              := scantokens (opt[nopt]arg);
          else:
            vlist[nfrom]arc[vlist[nfrom]arc.last]tag := 0;
          fi
        elseif match_option (opt[nopt], "width"):
          get_argument (opt[nopt], scantokens (opt[nopt]arg),
                        vlist[nfrom]arc[vlist[nfrom]arc.last]wd);
        elseif match_option (opt[nopt], "foreground"):
          get_argument (opt[nopt], scantokens (opt[nopt]arg),
                        vlist[nfrom]arc[vlist[nfrom]arc.last]fore);
        elseif match_option (opt[nopt], "background"):
          get_argument (opt[nopt], scantokens (opt[nopt]arg),
                        vlist[nfrom]arc[vlist[nfrom]arc.last]back);
        elseif match_option (opt[nopt], "rubout"):
          if known opt[nopt]arg:
             vlist[nfrom]arc[vlist[nfrom]arc.last]rub
               := scantokens (opt[nopt]arg);
          else:
             vlist[nfrom]arc[vlist[nfrom]arc.last]rub := 2;
          fi
        else:
          ignore_option (opt[nopt], opt[nopt]arg);
        fi
      endfor
      handle_line_style (vlist[nfrom]arc[vlist[nfrom]arc.last]sty, sty);
      vlist[nto]arc[vlist[nto]arc.last]lsr
        := vlist[nfrom]arc[vlist[nfrom]arc.last]lsr;
    fi
    nfrom := nto;
  endfor
enddef;
vardef handle_line_style (suffix sty) (expr name) =
  if valid_style name:
    sty := name;
  else:
    errhelp "feynmf: your linestyle is not recognizable, "
          & "check spelling and reprocess!";
    errmessage "feynmf: line style `" & name & "' not known, "
             & "replaced by `vanilla'";
    sty := "vanilla";
  fi
enddef;
vardef get_argument (expr opt, arg) (suffix variable) =
  if known arg:
    variable := arg;
  else:
    message "feynmf: option `" & opt & "' needs an argument.  Ignored.";
  fi
enddef;
vardef ignore_argument (expr opt, arg) =
  if known arg:
    message "feynmf: option `" & opt & "' doesn't take an argument.  "
          & "Argument `" & arg & "' ignored.";
  fi
enddef;
vardef ignore_option (expr opt, arg)=
  if known arg:
    message "feynmf: ignoring option " & opt & "=" & arg & ".";
  else:
    message "feynmf: ignoring option " & opt & ".";
  fi
enddef;
vardef vconnectn (expr linesty) (suffix v) (expr n) =
  vconnect (linesty, vmklist (v, n));
enddef;
vardef vpath@# (suffix from, to) =
  save nfrom, nto, origin, index, unknown_path;
  numeric nfrom, nto, origin, index;
  path unknown_path;
  if (known vloc from) and (known vloc to):
    nfrom := vlookup from;
    nto := vlookup to;
    vmatch_path (nfrom, nto, maybe_empty@#);
    if (unknown origin) or (unknown index):
      vmatch_path (nto, nfrom, maybe_empty@#);
    fi
  fi
  if (known origin) and (known index):
    vbuild_cut_arc (origin, index)
  else:
    unknown_path
  fi
enddef;
vardef maybe_empty@# =
  save _prefix;
  _prefix=137;
  if known _prefix@#:
    whatever
  else:
    @#
  fi
enddef;
vardef vmatch_path (expr nfrom, nto, t) =
  save i;
  for i = varcs (nfrom):
    if (vlist[nfrom]arc[i] = nto) and (known vlist[nfrom]arc[i]sty):
      if unknown t:
        origin := nfrom;
        index := i;
      else:
        if known vlist[nfrom]arc[i]tag:
          if vlist[nfrom]arc[i]tag = t:
            origin := nfrom;
            index := i;
          fi
        fi
      fi
    fi
  endfor
enddef;
vardef vcyclen (expr sty) (suffix v) (expr n) =
  for $ = 1 upto n - 1:
    vconnect (sty, v[$], v[$+1]);
  endfor
  vconnect (sty, v[n], v[1]);
enddef;
vardef vrcyclen (expr sty) (suffix v) (expr n) =
  vconnect (sty, v[1], v[n]);
  for $ = n downto 2:
    vconnect (sty, v[$], v[$-1]);
  endfor
enddef;
vardef vforce (expr z) (suffix v) =
  venter v;
  vlist[vlookup v]loc := z;
enddef;
vardef vshift (expr z) (text vl) =
  forsuffixes $=vl:
    if vexists $:
      vlist[vlookup $]loc := vlist[vlookup $]loc + z;
    fi
  endfor
enddef;
vardef vconstraint (expr z) (text vl) =
  save nfrom, nto;
  numeric nfrom, nto;
  forsuffixes to = vl:
    venter to;
    nto := vlookup to;
    if known nfrom:
      vlist[nfrom]constr.last := vlist[nfrom]constr.last + 1;
      vlist[nto]constr.last := vlist[nto]constr.last + 1;
      vlist[nfrom]constr[vlist[nfrom]constr.last] := nto;
      vlist[nto]constr[vlist[nto]constr.last] := nfrom;
      vlist[nto]loc = vlist[nfrom]loc + z;
    fi
    nfrom := nto;
  endfor
enddef;
vardef vpolygon (expr psty) (suffix v) (text vl) =
  save nopt, csty, nfrom, nfrom_, nto, i, n, j;
  numeric nopt, nfrom, nfrom_, nto, i, n, j;
  string csty;
  n := count (vl) + 1;
  plist.last := plist.last + 1;
  plist[plist.last]cnt := n;
  plist[plist.last]lbl := "";
  plist[plist.last]lbl.ang := whatever;
  plist[plist.last]lbl.dist := 0;
  csty := "phantom";
  getopt (opt, psty);
  for nopt = opt.first upto opt.last:
    if match_option (opt[nopt], "filled"):
      get_argument (opt[nopt], opt[nopt]arg, plist[plist.last]sty);
    elseif match_option (opt[nopt], "tension"):
      if known opt[nopt]arg:
        csty := csty & ",tension=" & opt[nopt]arg;
      else:
        message "feynmf: option `tension' needs an argument.  Ignored.";
      fi
    elseif match_option (opt[nopt], "label"):
      get_argument (opt[nopt], opt[nopt]arg, plist[plist.last]lbl);
    elseif match_option (opt[nopt], "label.angle"):
      get_argument (opt[nopt], scantokens (opt[nopt]arg),
                               plist[plist.last]lbl.ang);
    elseif match_option (opt[nopt], "label.dist"):
      get_argument (opt[nopt], scantokens (opt[nopt]arg),
                               plist[plist.last]lbl.dist);
    elseif match_option (opt[nopt], "pull"):
      get_argument (opt[nopt], scantokens (opt[nopt]arg),
                               plist[plist.last]pull);
    elseif match_option (opt[nopt], "cona"):
      get_argument (opt[nopt], opt[nopt]arg, plist[plist.last]cona);
    elseif match_option (opt[nopt], "conb"):
      get_argument (opt[nopt], opt[nopt]arg, plist[plist.last]conb);
    elseif match_option (opt[nopt], "smooth"):
      plist[plist.last]cona := "..";
      plist[plist.last]conb := "..";
      ignore_argument (opt[nopt], opt[nopt]arg);
    elseif match_option (opt[nopt], "foreground"):
      get_argument (opt[nopt], scantokens (opt[nopt]arg),
                               plist[plist.last]fore);
    elseif match_option (opt[nopt], "background"):
      get_argument (opt[nopt], scantokens (opt[nopt]arg),
                               plist[plist.last]back);
    elseif match_option (opt[nopt], "phantom"):
      plist[plist.last]sty := "phantom";
    elseif match_option (opt[nopt], "empty"):
      plist[plist.last]sty := "empty";
    elseif match_option (opt[nopt], "full"):
      plist[plist.last]sty := "full";
    elseif match_option (opt[nopt], "hatched"):
      plist[plist.last]sty := "hatched";
    elseif match_option (opt[nopt], "shaded"):
      plist[plist.last]sty := "shaded";
    else:
      ignore_option (opt[nopt], opt[nopt]arg);
  fi
  endfor
  canonicalize_filling plist[plist.last]sty;
  vconnect (csty, v, vl, v);
  i := 1;
  forsuffixes to = v, vl, v:
    nto := vlookup to;
    if known nfrom:
      if known nfrom_:
        vlist[nto]loc = vlist[nfrom]loc
          + (vlist[nfrom]loc - vlist[nfrom_]loc) rotated (360/n);
      fi
      vlist[nto]poly.last := vlist[nto]poly.last + 1;
      vlist[nto]poly[vlist[nto]poly.last] := plist.last;
      vlist[nto]poly[vlist[nto]poly.last]idx := i;
      plist[plist.last]vtx[i] := nto;
      i := i + 1;
      nfrom_ := nfrom;
    fi
    nfrom := nto;
  endfor
enddef;
vardef canonicalize_filling suffix f =
  if known f:
    if is_a_number (f):
      if scantokens f <= 1:
        f :=
        if scantokens f = 1:
          "full"
        elseif scantokens f > 0:
          "shaded"
        elseif scantokens f = 0:
          "empty"
        else:
          "hatched"
        fi;
      fi
    fi
  fi
enddef;
vardef vpolygonn (expr sty) (suffix v) (expr n) =
  vpolygon (sty, v[1], for $=2 upto n-1: v[$], endfor v[n]);
enddef;
vardef vrpolygonn (expr sty) (suffix v) (expr n) =
  vpolygon (sty, v[n], for $=n-1 downto 2: v[$], endfor v[1]);
enddef;
vardef vlabel (expr s) (suffix v) =
  venter v;
  vlist[vlookup v]lbl := s;
enddef;
vardef vvertex (expr vtxsty) (text vl) =
  save nopt, sty, arg;
  numeric nopt, arg;
  string sty;
  getopt (opt, vtxsty);
  forsuffixes v = vl:
    venter v;
    n := vlookup v;
    for nopt = opt.first upto opt.last:
      handle_vertex_option (vlist[n], opt[nopt], opt[nopt]arg);
    endfor
  endfor
enddef;
vardef handle_vertex_option (suffix v) (expr opt, arg) =
  if match_option (opt, "label"):
    get_argument (opt, arg, v.lbl);
  elseif match_option (opt, "label.angle"):
    get_argument (opt, scantokens (arg), v.lbl.ang);
  elseif match_option (opt, "label.dist"):
    get_argument (opt, scantokens (arg), v.lbl.dist);
  elseif match_option (opt, "decoration.shape"):
    if known arg:
      make_decor_shape (v.decor.shape, arg);
    else:
      message "feynmf: option `decor.shape' needs an argument.  Ignored.";
    fi
  elseif match_option (opt, "decoration.filled"):
    get_argument (opt, arg, v.decor.sty);
    canonicalize_filling v.decor.sty;
  elseif match_option (opt, "decoration.size"):
    get_argument (opt, scantokens (arg), v.decor.size);
  elseif match_option (opt, "decoration.angle"):
    get_argument (opt, scantokens (arg), v.decor.ang);
  elseif match_option (opt, "foreground"):
    get_argument (opt, scantokens (arg), v.fore);
  elseif match_option (opt, "background"):
    get_argument (opt, scantokens (arg), v.back);
  else:
    ignore_option (opt, arg);
  fi
enddef;
vardef make_decor_shape (suffix p) (expr n) =
  if match_prefix (n, "circle"): p := fullcircle;
  elseif match_prefix (n, "square"):
    p := unitsquare shifted -(.5,.5);
  elseif match_prefix (n, "cross"): p := polycross 4;
  elseif match_prefix (n, "triangle"): p := polygon 3;
  elseif match_prefix (n, "triagon"): p := polygon 3;
  elseif match_prefix (n, "diamond"): p := polygon 4;
  elseif match_prefix (n, "tetragon"): p := polygon 4;
  elseif match_prefix (n, "pentagon"): p := polygon 5;
  elseif match_prefix (n, "hexagon"): p := polygon 6;
  elseif match_prefix (n, "triagram"): p := polygram 3;
  elseif match_prefix (n, "tetragram"): p := polygram 4;
  elseif match_prefix (n, "pentagram"): p := polygram 5;
  elseif match_prefix (n, "hexagram"): p := polygram 6;
  elseif match_prefix (n, "triacross"): p := polycross 3;
  elseif match_prefix (n, "tetracross"): p := polycross 4;
  elseif match_prefix (n, "pentacross"): p := polycross 5;
  elseif match_prefix (n, "hexacross"): p := polycross 6;
  else:
    if feynmfwizard:
      p := scantokens(n);
    else:
      message "feynmf: invalid argument `" & n
            & "' to option `decor.shape'.  Ignored.";
    fi
  fi
enddef;
vardef is_a_number expr s =
  save n;
  if known s:
    (true
     for n = 1 upto length s:
       and ((substring (n-1,n) of s = ".")
            or (substring (n-1,n) of s = "-")
            or isdigit substring (n-1,n) of s)
      endfor)
  else:
    false
  fi
enddef;
vardef vvertexn (expr vtxsty) (suffix v) (expr n) =
  vvertex (vtxsty, vmklist (v, n));
enddef;
vardef vblob (expr bd) (text vl)=
  forsuffixes $=vl:
    if not vexists $: venter $; fi
    vlist[vlookup $]decor.shape := fullcircle;
    vlist[vlookup $]decor.size := bd;
    vlist[vlookup $]decor.sty := "shaded";
 endfor
enddef;
vardef vdot (text vl)=
  forsuffixes $=vl:
    if not vexists $: venter $; fi
    vlist[vlookup $]decor.shape := fullcircle;
    vlist[vlookup $]decor.size := dot_size;
    vlist[vlookup $]decor.sty := "full";
 endfor
enddef;
vardef vdotn (suffix v) (expr n) =
  vdot (vmklist (v, n));
enddef;
vardef vblobn (suffix v) (expr n) =
  vblob (vmklist (v, n));
enddef;
vardef curved_left_gallery = .9[se,sw] .. .5[sw,nw] .. .1[nw,ne] enddef;
vardef curved_right_gallery = .9[sw,se] .. .5[se,ne] .. .1[ne,nw] enddef;
vardef curved_bottom_gallery = .9[nw,sw] .. .5[sw,se] .. .1[se,ne] enddef;
vardef curved_top_gallery = .9[sw,nw] .. .5[nw,ne] .. .1[ne,se] enddef;
vardef curved_surround_gallery =
  superellipse (.5[se,ne], .5[ne,nw], .5[nw,sw], .5[sw,se], .75)
enddef;
vardef straight_left_gallery = sw -- nw enddef;
vardef straight_right_gallery = se -- ne enddef;
vardef straight_bottom_gallery = sw -- se enddef;
vardef straight_top_gallery = nw -- ne enddef;
vardef straight_surround_gallery =
   .5[se,ne] -- ne -- .5[ne,nw] -- nw
     -- .5[nw,sw] -- sw -- .5[sw,se] -- se -- cycle
enddef;
vardef curved_galleries =
  vardef left_gallery = curved_left_gallery enddef;
  vardef right_gallery = curved_right_gallery enddef;
  vardef bottom_gallery = curved_bottom_gallery enddef;
  vardef top_gallery = curved_top_gallery enddef;
  vardef surround_gallery = curved_surround_gallery enddef;
enddef;
vardef straight_galleries =
  vardef left_gallery = straight_left_gallery enddef;
  vardef right_gallery = straight_right_gallery enddef;
  vardef bottom_gallery = straight_bottom_gallery enddef;
  vardef top_gallery = straight_top_gallery enddef;
  vardef surround_gallery = straight_surround_gallery enddef;
enddef;
vardef vleft (text vl) = vdistribute (left_gallery, vl) enddef;
vardef vright (text vl) = vdistribute (right_gallery, vl) enddef;
vardef vbottom (text vl) = vdistribute (bottom_gallery, vl) enddef;
vardef vtop (text vl) = vdistribute (top_gallery, vl) enddef;
vardef vsurround (text vl) = vdistribute (surround_gallery, vl) enddef;
curved_galleries;
vardef vdistribute (expr p) (text vl) =
  save numv, len, off;
  numeric numv, len, off;
  numv = count (vl);
  if cycle p: numv := numv + 1; fi
  len := length (p);
  if numv = 1:
    vforce (point len/2 of p, vl);
  else:
    off := 0;
    forsuffixes $ = vl:
      vforce (point off of p, $);
      off := off + len/(numv-1);
    endfor
  fi
enddef;
def vmklist (suffix v) (expr n) =
  for $ = 1 upto n-1: v[$], endfor v[n]
enddef;
vardef vleftn (suffix v) (expr n) =
  vleft (vmklist (v, n));
enddef;
vardef vrightn (suffix v) (expr n) =
  vright (vmklist (v, n));
enddef;
vardef vbottomn (suffix v) (expr n) =
  vbottom (vmklist (v, n));
enddef;
vardef vtopn (suffix v) (expr n) =
  vtop (vmklist (v, n));
enddef;
vardef vsurroundn (suffix v) (expr n) =
  vsurround (vmklist (v, n));
enddef;
vardef vfreeze =
  for i = vertices:
    if unknown vlist[i]loc:
      origin = origin
      for j = varcs (i):
        + vlist[i]arc[j]tns * (vlist[i]loc - vlist[vlist[i]arc[j]]loc)
      endfor
      for j = vconstr (i):
        if i < vlist[i]constr[j]:
          + lambda (i, vlist[i]constr[j])
        elseif i > vlist[i]constr[j]:
          - lambda (vlist[i]constr[j], i)
        fi
      endfor
      for j = vpoly (i):
        + lambdapoly (vlist[i]poly[j], plist[vlist[i]poly[j]]cnt,
                      vlist[i]poly[j]idx)
      endfor;
    fi
  endfor
  if vtracing: vdump; fi
enddef;
vardef lambda (expr i, j) =
  (if known (xpart(vlist[i]loc - vlist[j]loc)):
     lambdax[i][j]
   else:
     0
   fi,
   if known (ypart(vlist[i]loc - vlist[j]loc)):
     lambday[i][j]
   else:
     0
   fi)
enddef;
vardef lambdapoly (expr p, n, i) =
  origin
  if i = 1:
    + lambdap[p][2] rotated (-360/n)
  elseif i = 2:
    - lambdap[p][1]
  fi
  if i = n - 1:
    - lambdap[p][n] rotated (-360/n)
  elseif i = n:
    + lambdap[p][n-1]
  fi
  if (i > 1) and (i < n):
    + lambdap[p][i-1] + lambdap[p][i+1] rotated (-360/n)
    - lambdap[p][i] - lambdap[p][i] rotated (-360/n)
  fi
enddef;
vardef idraw (expr linesty, p) =
  save nopt, sty, lbl, wd, rub;
  numeric nopt, lbl.dist, wd, rub;
  save  fore, back;
  color fore, back;
  string sty, lbl, lbl.side;
  getopt (opt, linesty);
  sty := opt[opt.first];
  if known opt[opt.first]arg:
    message "feynmf: line styles don't take arguments.  "
             & "Argument `" & opt[opt.first]arg & "' ignored.";
  fi
  opt.first := opt.first + 1;
  lbl := "";
  lbl.side := "";
  lbl.dist := 3thick;
  for nopt = opt.first upto opt.last:
    if match_option (opt[nopt], "label"):
      get_argument (opt[nopt], opt[nopt]arg, lbl);
    elseif match_option (opt[nopt], "label.side"):
      get_argument (opt[nopt], opt[nopt]arg, lbl.side);
    elseif match_option (opt[nopt], "label.dist"):
      get_argument (opt[nopt], scantokens (opt[nopt]arg), lbl.dist);
    elseif match_option (opt[nopt], "width"):
      get_argument (opt[nopt], scantokens (opt[nopt]arg), wd);
    elseif match_option (opt[nopt], "foreground"):
      get_argument (opt[nopt], scantokens (opt[nopt]arg), fore);
    elseif match_option (opt[nopt], "background"):
      get_argument (opt[nopt], scantokens (opt[nopt]arg), back);
    elseif match_option (opt[nopt], "rubout"):
      if known opt[nopt]arg:
         rub := scantokens (opt[nopt]arg);
      else:
         rub := 2;
      fi
    else:
      ignore_option (opt[nopt], opt[nopt]arg);
    fi
  endfor
  handle_line_style (sty, sty);
  begingroup
    if known fore:
      save foreground;
      color foreground;
      foreground = fore;
    fi
    if known back:
      save background;
      color background;
      foreground = back;
    fi
    vdraw_arc_rubout (rub, sty, wd, p, lbl);
  endgroup;
enddef;
vardef ivertex (expr vtxsty, pos) =
  save nopt, v;
  numeric nopt, v.lbl.ang, v.lbl.dist,
    v.decor.size, v.decor.ang;
  pair v.loc;
  string v.lbl, v.decor.sty;
  path v.decor.shape;
  v.loc := pos;
  v.lbl := "";
  v.lbl.ang := whatever;
  v.lbl.dist := 3thick;
  v.decor.size := decor_size;
  getopt (opt, vtxsty);
  for nopt = opt.first upto opt.last:
    handle_vertex_option (v, opt[nopt], opt[nopt]arg);
  endfor
  vdraw_label (v.loc, v.lbl);
  vdraw_vertex v;
enddef;
vardef vdraw =
  if not feynmfwizard:
    vcheck_typos;
  fi
  for i = vertices:
    if not known vlist[i]loc:
      errhelp "Your graph specification was not complete (probably a "
            & "lone vertex).            Check logic and reprocess!";
      errmessage "feynmf: vertex `" & vlist[i]name & "' not determined, "
               & "replaced by `(0,0)'.";
      vlist[i]loc := origin;
    fi
    if unknown vlist[i]decor.size:
      vlist[i]decor.size = decor_size;
    fi
  endfor
  for i = vertices:
    for j = varcs (i):
      if known vlist[i]arc[j]sty:
        if vlist[i]arc[j]sty <> "":
          if unknown vlist[i]arc[j]rub:
            begingroup
              if known vlist[i]arc[j]fore:
                save foreground; color foreground;
                foreground = vlist[i]arc[j]fore;
              fi
              if known vlist[i]arc[j]back:
                save background; color background;
                background = vlist[i]arc[j]back;
              fi
              vdraw_arc (vlist[i]arc[j]sty, vlist[i]arc[j]wd,
                         vbuild_cut_arc (i, j), vlist[i]arc[j]lbl);
              vlist[i]arc[j]sty := "";
            endgroup;
          fi
        fi
      fi
    endfor;
  endfor
  for i = polygons:
    vdraw_label (pcenter plist[i], plist[i]lbl);
    vdraw_polygon plist[i];
  endfor
  for i = vertices:
    vdraw_label (vlist[i]loc, vlist[i]lbl);
    vdraw_vertex vlist[i];
  endfor
  for i = vertices:
    for j = varcs (i):
      if known vlist[i]arc[j]sty:
        if vlist[i]arc[j]sty <> "":
          if known vlist[i]arc[j]rub:
            begingroup
              if known vlist[i]arc[j]fore:
                save foreground; color foreground;
                foreground = vlist[i]arc[j]fore;
              fi
              if known vlist[i]arc[j]back:
                save background; color background;
                background = vlist[i]arc[j]back;
              fi
              vdraw_arc_rubout (vlist[i]arc[j]rub,
                                vlist[i]arc[j]sty,
                                vlist[i]arc[j]wd,
                                vbuild_cut_arc (i, j),
                                vlist[i]arc[j]lbl);
              vlist[i]arc[j]sty := "";
            endgroup;
          fi
        fi
      fi
    endfor;
  endfor
enddef;
vardef vcheck_typos =
  save j, err;
  boolean wrn;
  wrn := false;
  for i = vertices:
    save connections;
    connections = vlist[i]arc.last - vlist[i]arc.first + 1;
    if connections < 1:
      if unknown vlist[i]loc:
        message "feynmf: warning: disconnected and unspecified vertex `"
              & substring (2,infinity) of vlist[i]name
              & "'.";
        wrn := true;
      fi
    elseif connections = 1:
      j := vlist[i]arc[vlist[i]arc.last];
      if j < i:
        if vlist[i]loc = vlist[j]loc:
          message "feynmf: warning: dangling vertex `"
                & substring (2,infinity) of vlist[i]name
                & "' colliding with `"
                & substring (2,infinity) of vlist[j]name
                & "'.";
        wrn := true;
        fi
      fi
    fi
  endfor
  if wrn:
    message "feynmf: Have you seen the warning messages above?";
    message "        They are usually caused by misspelling a vertex'";
    message "        name and can trigger errors further below!";
    message "        Fix the typos and run LaTeX and Metafont again.";
  fi
enddef;
vardef vbuild_arc (expr lsr, from, to) =
  if unknown lsr:
    from -- to
  else:
    if lsr = 0:
      from -- to
    else:
      from
        .. (1-lsr)/2 *(to rotatedabout (.5[from,to], 90))
           + (1+lsr)/2 * (to rotatedabout (.5[from,to], -90))
        .. to
    fi
  fi
enddef;
vardef vbuild_cut_arc (expr origin, index) =
  cut_decors (vlist[origin],
              if vlist[origin]arc[index] <> origin:
                vbuild_arc (vlist[origin]arc[index]lsr,
                            vlist[origin]loc,
                            vlist[vlist[origin]arc[index]]loc)
              else:
                vbuild_tadpole (origin, index)
              fi,
              vlist[vlist[origin]arc[index]])
enddef;
vardef vbuild_tadpole (expr origin, index) =
  save j, n, nn, nnn, aidx, aang, agap, bgap, ngap, distsum;
  n := 0;
  distsum := 0;
  for j = varcs (origin):
    if vlist[origin]arc[j] <> origin:
      ang := angle direction 0 of
        vbuild_arc (vlist[origin]arc[j]lsr,
                    vlist[origin]loc,
                    vlist[vlist[origin]arc[j]]loc);
      n := n + 1;
      distsum := distsum
        + abs (vlist[vlist[origin]arc[j]]loc - vlist[origin]loc);
      aang[n] := 360;
      for nn = 1 upto n:
        if ang < aang[nn]:
          for nnn = n - 1 downto nn:
            aang[nnn+1] := aang[nnn];
            aidx[nnn+1] := aidx[nnn];
          endfor
          aang[nn] := ang;
          aidx[nn] := n;
        fi
        exitif known aidx[n];
      endfor
    fi
  endfor
  aidx[n+1] := aidx[1];
  aang[n+1] := aang[1] + 360;
  for nn = 1 upto n:
    agap[nn] = aang[nn+1] - aang[nn];
  endfor
  if known vlist[origin]arc[index]lsr:
    ngap := n;
    for nn = 1 upto n:
      if (aang[nn] < vlist[origin]arc[index]lsr)
         and (vlist[origin]arc[index]lsr < aang[nn+1]):
        ngap := nn;
      fi
    endfor
  else:
    bgap := 0;
    for nn = 1 upto n:
      if agap[nn] > bgap:
        bgap :=  agap[nn];
        ngap := nn;
      fi
    endfor
  fi
  if vtracing: adump (n + 1); fi
  vlist[origin]loc{dir(aang[ngap]+agap[ngap]/4)}
    ... vlist[origin]loc + 2/3 * distsum/n
                            / vlist[origin]arc[index]tns
                            * dir(aang[ngap]+agap[ngap]/2)
    ... {-dir(aang[ngap+1]-agap[ngap]/4)}vlist[origin]loc
enddef;
vardef adump expr n =
  save i;
  for i = 1 upto n:
    message "aidx[" & decimal_ (i) & "]=" & decimal_ (aidx[i])
            & ", aang[" & decimal_ (i) & "]=" & decimal_ (aang[i])
            & ", agap[" & decimal_ (i) & "]=" & decimal_ (agap[i]);
  endfor
enddef;
vardef vdraw_arc (expr sty, wd, arc) (suffix lbl) =
  if known wd:
    save currentpen;
    pen currentpen;
    pickup pencircle scaled wd;
  fi
  scantokens ("draw_" & sty) (arc);
  vdraw_arc_label (arc, lbl);
enddef;
let plain_draw = draw;
vardef vdraw_arc_rubout (expr rub, sty, wd, arc) (suffix lbl) =
  if known rub:
    begingroup
      def draw expr p =
        save oldpen; pen oldpen;
        oldpen := currentpen;
        save currentpen; pen currentpen;
        pickup oldpen scaled rub;
        erase plain_draw (subpath (.1,.9)*length(p) of p)
      enddef;
      vdraw_arc (sty, wd, arc, lbl);
      let draw = plain_draw;
    endgroup;
  fi
  vdraw_arc (sty, wd, arc, lbl);
enddef;
vardef vbuild_polygon suffix p =
  if known p.pull:
    save c; pair c;
    c := pcenter p;
    for i = 1 upto (p.cnt - 1):
      vbuild_polygon_section (p, i, i+1)
    endfor
    vbuild_polygon_section (p, p.cnt, 1)
  else:
    for i = 1 upto p.cnt:
      vlist[p.vtx[i]]loc
        if known p.cona: scantokens (p.cona) else: -- fi
    endfor
  fi
  cycle
enddef;
def vbuild_polygon_section (suffix p) (expr from, to) =
  vlist[p.vtx[from]]loc
    if known p.cona: scantokens (p.cona) else: -- fi
  p.pull[c,.5[vlist[p.vtx[from]]loc,vlist[p.vtx[to]]loc]]
    if known p.conb: scantokens (p.conb) else: -- fi
enddef;
vardef pcenter suffix p =
 (origin for i = 1 upto p.cnt:
           + vlist[p.vtx[i]]loc
         endfor) / p.cnt
enddef;
vardef vdraw_arc_label (expr arc) (suffix lbl) =
  if lbl <> "":
    save _a, _z, _zz, _r;
    numeric _a;
    pair _z, _zz, _r;
    _z := point .5 length (arc) of arc;
    if lbl.dist = 0:
      LaTeX_text (_z, whatever, lbl);
    else:
      _r := direction .5 length (arc) of arc rotated - 90;
      if lbl.side = "left":
        _a := angle (-_r);
      elseif lbl.side = "right":
        _a := angle (_r);
      else:
        _zz = _z - .5[point 0 of arc, point infinity of arc];
        if ((_zz if length (_zz) >  0: / length (_zz) fi)
            dotprod _r) >= 0:
          _a := angle (_r);
        else:
          _a := angle (-_r);
        fi
      fi
      LaTeX_text (_z + lbl.dist * dir _a, _a, lbl);
    fi
  fi
enddef;
vardef vdraw_label (expr loc) (suffix lbl) =
  if lbl <> "":
    save a;
    numeric a;
    if lbl.dist = 0:
      LaTeX_text (loc, whatever, lbl);
    else:
      if unknown lbl.ang:
        if loc = (.5w,.5h):
          a := 0;
        else:
          a := angle (loc - (.5w,.5h));
        fi
      else:
        a := lbl.ang;
      fi
      LaTeX_text (loc + lbl.dist * dir a, a, lbl);
    fi
  fi
enddef;
vardef vdraw_vertex suffix v =
  if known v.decor.shape:
    if known v.fore:
      save foreground;
      color foreground;
      foreground = v.fore;
    fi
    if known v.back:
      save background;
      color background;
      background = v.back;
    fi
    if known v.decor.sty:
      if v.decor.sty = "empty":
        emptydraw (
      elseif v.decor.sty = "full":
        cfilldraw (
      elseif is_a_number v.decor.sty:
        drawhalftone (1 - scantokens v.decor.sty / 100,
      else:
        drawtile (scantokens v.decor.sty,
      fi
    else:
      cfilldraw (
    fi
        v.decor.shape
          if known v.decor.ang: rotated v.decor.ang fi
          scaled v.decor.size shifted v.loc);
  fi
enddef;
vardef polygon expr n =
  if n > 2:
    for i = 1 upto n:
      (.5up rotated (360i/n)) --
    endfor
    cycle
  else:
    fullcircle
  fi
enddef;
vardef polygram expr n =
  if n > 2:
    for i = 1 upto n:
      (.5up rotated (360i/n)) --
      (.2up rotated (360(i+.5)/n)) --
    endfor
    cycle
  else:
    fullcircle
  fi
enddef;
vardef polycross expr n =
  save i;
  for i = 1 upto n:
    origin -- .5 dir (360(i-.5)/n) --
  endfor
  cycle
enddef;
vardef vdraw_polygon suffix p =
  if known p.fore:
    save foreground;
    color foreground;
    foreground = p.fore;
  fi
  if known p.back:
    save background;
    color background;
    background = p.back;
  fi
  if known p.sty:
    if p.sty = "phantom":
    elseif p.sty = "empty":
      emptydraw (vbuild_polygon p);
    elseif p.sty = "full":
      cfilldraw (vbuild_polygon p);
    elseif is_a_number p.sty:
      drawhalftone (1 - scantokens p.sty / 100, vbuild_polygon p);
    else:
      drawtile (scantokens p.sty, vbuild_polygon p);
    fi
  else:
    cdraw (vbuild_polygon p);
  fi
enddef;
vardef LaTeX expr text =
  if LaTeX_file = "":
    LaTeX_file := jobname & ".t" & decimal charcode;
    write ("% " & LaTeX_file & " -- generated from " & jobname & ".mp")
      to LaTeX_file;
  fi
  write (text & "%") to LaTeX_file
enddef;
vardef LaTeX_text (expr z, a, txt) =
  LaTeX "\fmfL(" & (decimal (xpart z/LaTeX_unitlength)) & ","
      & (decimal (ypart z/LaTeX_unitlength)) & ","
      & (voctant a) & "){" & txt & "}";
enddef;
vardef voctant expr a =
  if known a:
    voctant_list[floor (a/45 + .5)]
  else:
    "c"
  fi
enddef;
string voctant_list[];
voctant_list[-4] := "r";
voctant_list[-3] := "rt";
voctant_list[-2] := "t";
voctant_list[-1] := "lt";
voctant_list[0] := "l";
voctant_list[1] := "lb";
voctant_list[2] := "b";
voctant_list[3] := "rb";
voctant_list[4] := "r";
vardef vdump =
  message ">>>>> Vertices and arcs for diagram #" & decimal charcode
        & " of " & jobname & ".mf:";
  for i = vertices:
    message "> " & vlist[i]name & "=" & decimal_pair (vlist[i]loc)
          & ": #lines="
          & decimal (vlist[i]arc.last - vlist[i]arc.first + 1)
          if vlist[i]lbl <> "":
            & ", lbl=" & vlist[i]lbl
            & ", l.angle=" & decimal_ (vlist[i]lbl.ang)
            & ", l.dist=" & decimal_ (vlist[i]lbl.dist)
          fi
          & ".";
  endfor
  for i = vertices:
    for j = varcs (i):
      if known vlist[i]arc[j]sty:
        message "> " & vlist[i]name & "*" & vlist[vlist[i]arc[j]]name
                & ": " & vlist[i]arc[j]sty
                & ", tns=" & decimal_ (vlist[i]arc[j]tns)
                & ", lsr=" & decimal_ (vlist[i]arc[j]lsr)
                & ", wd=" & decimal_ (vlist[i]arc[j]wd)
                & ", rub=" & decimal_ (vlist[i]arc[j]rub)
                if vlist[i]arc[j]lbl <> "":
                  & ", lbl=" & vlist[i]arc[j]lbl
                  & ", l.side=" & vlist[i]arc[j]lbl.side
                  & ", l.dist=" & decimal_ (vlist[i]arc[j]lbl.dist)
                fi
                & ".";
      fi
    endfor
    for j = vconstr (i):
      if i < vlist[i]constr[j]:
        save z;
        pair z;
        z = vlist[vlist[i]constr[j]]loc - vlist[i]loc;
        message "> " & vlist[i]name & "&"
                & vlist[vlist[i]constr[j]]name
                & ": " & decimal_pair (z);
      fi
    endfor;
  endfor
enddef;
vardef decimal_ (text n) =
  if known n: decimal n else: "?" fi
enddef;
vardef decimal_pair (text z) =
  "(" & decimal_ (xpart z) & "," & decimal_ (ypart z) & ")"
enddef;
def show_diagram_ expr frame =
  if (screen_cols < w + 2 xpart frame) or (screen_rows < h + 2 ypart frame):
    screen_cols := w + 2 xpart frame;
    screen_rows := h + 2 ypart frame;
    openwindow currentwindow
      from origin to (screen_rows, screen_cols)
      at (- xpart frame, h + ypart frame);
  fi
  showit_;
  if showstopping > 0:
    stop "This is diagram #" & decimal charcode
       & ".  Hit return to continue...";
  fi
enddef;
def show_diagram =
  def show_diagram =
    display blankpicture inwindow currentwindow;
    show_diagram_
  enddef;
  show_diagram_
enddef;
def show_all_diagrams expr frame =
  def showit = show_diagram frame enddef;
  displaying:=1;
enddef;
endinput;
\endinput
%%
%% End of file `feynmp.mp'.
