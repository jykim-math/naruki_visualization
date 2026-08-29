const state={lambda:2,mu:3,nu:5,rho:.25,points:[[],[],[],[]],views:Array.from({length:4},()=>({yaw:-.62,pitch:.48,zoom:1})),selected:null,lineModels:new Map()};
const COLORS={AB:"#ed775f",AC:"#668bb3",BC:"#789467"};
const LINE_COLORS=Array.from({length:27},(_,i)=>"hsl("+((i*137.508)%360)+" 48% 61%)");
const LINES=[
["a1","B1","C1"],["a2","B3","C3"],["a3","B2","C2"],["a4","B3","C2"],["a5","B2","C3"],["a6","A2","B1"],["a7","A2","C1"],["a8","A3","C1"],["a9","A3","B1"],
["b1","A1","C1"],["b2","A3","C3"],["b3","A2","C2"],["b4","A2","C3"],["b5","A3","C2"],["b6","B2","C1"],["b7","A1","B2"],["b8","A1","B3"],["b9","B3","C1"],
["c1","A1","B1"],["c2","A3","B3"],["c3","A2","B2"],["c4","A3","B2"],["c5","A2","B3"],["c6","A1","C2"],["c7","B1","C2"],["c8","B1","C3"],["c9","A1","C3"]
].map(([label,from,to],i)=>({label,from,to,family:[from[0],to[0]].sort().join(""),color:LINE_COLORS[i]}));
const endpoints=new Map(LINES.map(x=>[[x.from,x.to].sort().join("-"),x]));
const canvases=[...document.querySelectorAll("[data-chart]")],rhoInput=document.querySelector("#rho");
let timer;
const rhoZero=()=>Math.abs(state.rho)<1e-9;

function homogeneous([x0,x1,x2,x3],rho=state.rho){
 const {lambda:l,mu:m,nu:n}=state,q=l*m*n*rho-1,rm1=rho-1;
 const inner=l*x0*x0+m*x1*x1+n*x2*x2+(m*n+1)*x1*x2+(l*n+1)*x0*x2+(l*m+1)*x0*x1
  +rm1*rm1*q*q*x3*x3-rm1*q*x3*((l+1)*x0+(m+1)*x1+(n+1)*x2);
 return rho*x3*inner+x0*x1*x2;
}
function signedTerm(c,body){return (c<0?" − ":" + ")+Math.abs(c).toFixed(3)+body}
function updateEquation(){
 const {lambda:l,mu:m,nu:n,rho:r}=state,f=x=>Number(x).toFixed(3),square=(r-1)*(r-1)*(l*m*n*r-1)*(l*m*n*r-1),mixed=(r-1)*(l*m*n*r-1);
 document.querySelector("#surface-equation").textContent=
  f(r)+"x₃ [ "+f(l)+"x₀² + "+f(m)+"x₁² + "+f(n)+"x₂² + "+f(square)+"x₃² + "+
  f(m*n+1)+"x₁x₂ + "+f(l*n+1)+"x₀x₂ + "+f(l*m+1)+"x₀x₁"+
  signedTerm(-mixed,"x₃(("+f(l+1)+")x₀ + ("+f(m+1)+")x₁ + ("+f(n+1)+")x₂)")+" ] + x₀x₁x₂ = 0";
}
function chartValue(chart,u,v,w){const x=[],free=[0,1,2,3].filter(i=>i!==chart);x[chart]=1;x[free[0]]=u;x[free[1]]=v;x[free[2]]=w;return homogeneous(x)}
function gradient(chart,x,y,z){const h=.012;return[(chartValue(chart,x+h,y,z)-chartValue(chart,x-h,y,z))/(2*h),(chartValue(chart,x,y+h,z)-chartValue(chart,x,y-h,z))/(2*h),(chartValue(chart,x,y,z+h)-chartValue(chart,x,y,z-h))/(2*h)]}
function sampleChart(chart){
 const animated=typeof playing!=="undefined"&&playing.size>0,e=2.4,g=innerWidth<700?(animated?25:43):(animated?31:51),steps=animated?38:62,out=[],axes=[[0,1,2],[0,2,1],[1,2,0]];
 for(const ax of axes)for(let i=0;i<g;i++){const u=-e+2*e*i/(g-1);for(let j=0;j<g;j++){const v=-e+2*e*j/(g-1);let p=[0,0,0],old=-e;p[ax[0]]=u;p[ax[1]]=v;p[ax[2]]=old;let fo=chartValue(chart,...p);for(let k=1;k<=steps;k++){const t=-e+2*e*k/steps;p=[0,0,0];p[ax[0]]=u;p[ax[1]]=v;p[ax[2]]=t;const f=chartValue(chart,...p);if(f===0||fo*f<0){let lo=old,hi=t,flo=fo;for(let q=0;q<8;q++){const mid=(lo+hi)/2,test=[0,0,0];test[ax[0]]=u;test[ax[1]]=v;test[ax[2]]=mid;const fm=chartValue(chart,...test);if(flo*fm<=0)hi=mid;else{lo=mid;flo=fm}}const root=(lo+hi)/2,point=[0,0,0];point[ax[0]]=u;point[ax[1]]=v;point[ax[2]]=root;const n=gradient(chart,...point),nl=Math.hypot(...n)||1;out.push({p:point,n:n.map(a=>a/nl)});break}old=t;fo=f}}}
 return out;
}
function sample(){
 timer=null;
 document.querySelector("#surface-status").textContent="sampling four charts…";
 state.points=canvases.map((_,chart)=>sampleChart(chart));state.lineModels=computeLineModels();
 const total=state.points.reduce((n,p)=>n+p.length,0);document.querySelector("#sample-count").textContent=total.toLocaleString();document.querySelector("#viewing-state").textContent=rhoZero()?"boundary fiber · four charts":"four real affine loci";document.querySelector("#surface-status").textContent=rhoZero()?"boundary fiber":"four charts complete";drawAll();
}
function rotate(p,view){const[x,y,z]=p,cy=Math.cos(view.yaw),sy=Math.sin(view.yaw),cp=Math.cos(view.pitch),sp=Math.sin(view.pitch),xx=cy*x+sy*z,zz=-sy*x+cy*z;return[xx,cp*y-sp*zz,sp*y+cp*zz]}
function solve4(a,b){
 const m=a.map((r,i)=>[...r,b[i]]);
 for(let c=0;c<4;c++){let p=c;for(let i=c+1;i<4;i++)if(Math.abs(m[i][c])>Math.abs(m[p][c]))p=i;if(Math.abs(m[p][c])<1e-12)return null;[m[c],m[p]]=[m[p],m[c]];for(let i=c+1;i<4;i++){const q=m[i][c]/m[c][c];for(let j=c;j<5;j++)m[i][j]-=q*m[c][j]}}
 const x=[0,0,0,0];for(let i=3;i>=0;i--){x[i]=(m[i][4]-m[i].slice(i+1,4).reduce((s,v,j)=>s+v*x[i+1+j],0))/m[i][i]}return x;
}
function lineFrame(line,u){
 const piv={AB:[0,1],AC:[0,2],BC:[1,2]}[line.family],free=[0,1,2,3].filter(i=>!piv.includes(i)),p=[0,0,0,0],q=[0,0,0,0];p[piv[0]]=1;q[piv[1]]=1;p[free[0]]=u[0];p[free[1]]=u[1];q[free[0]]=u[2];q[free[1]]=u[3];return{p,q};
}
function lineResidual(line,u,rho){
 const {p,q}=lineFrame(line,u),f=t=>homogeneous(p.map((v,i)=>v+t*q[i]),rho),f0=f(0),f1=f(1),fm=f(-1),f2=f(2),c2=(f1+fm)*.5-f0,odd=(f1-fm)*.5,c3=(f2-f0-4*c2-2*odd)/6;return[f0,odd-c3,c2,c3];
}
function continueLine(line){
 const pa=boundaryPoint(line.from),pb=boundaryPoint(line.to);let u=[0,pa[3],0,pb[3]];if(rhoZero())return lineFrame(line,u);
 const target=state.rho,stages=[];for(let rho=Math.min(1e-6,target);rho<target;rho=Math.min(target,rho*1.2))stages.push(rho);if(stages.at(-1)!==target)stages.push(target);
 for(const rho of stages)for(let k=0;k<30;k++){const f=lineResidual(line,u,rho),norm=Math.max(...f.map(Math.abs));if(norm<1e-9)break;const h=1e-5,j=f.map(()=>[0,0,0,0]);for(let c=0;c<4;c++){const v=[...u];v[c]+=h;const g=lineResidual(line,v,rho);for(let r=0;r<4;r++)j[r][c]=(g[r]-f[r])/h}const d=solve4(j,f.map(x=>-x));if(!d||d.some(x=>!Number.isFinite(x)))return lineFrame(line,u);for(let i=0;i<4;i++)u[i]+=d[i]}
 return lineFrame(line,u);
}
function computeLineModels(){return new Map(LINES.map(line=>[line.label,continueLine(line)]))}
function boundaryPoint(label){
 const axis={A:0,B:1,C:2}[label[0]],index=+label[1],x=[0,0,0,0];
 x[axis]=1;
 if(index===2)x[3]=1;
 if(index===3)x[3]={A:state.lambda,B:state.mu,C:state.nu}[label[0]];
 return x;
}
function boundaryLineInChart(line,chart){
 const model=state.lineModels.get(line.label)||{p:boundaryPoint(line.from),q:boundaryPoint(line.to)},a=model.p,b=model.q,free=[0,1,2,3].filter(i=>i!==chart),e=2.4,samples=[];
 // A projective real line is a circle modulo antipodes. Angle coordinates also cover its point at infinity.
 for(let i=0;i<=360;i++){
  const t=Math.PI*i/360,c=Math.cos(t),s=Math.sin(t),h=a.map((v,j)=>c*v+s*b[j]),d=h[chart];
  if(Math.abs(d)<1e-7){samples.push(null);continue}
  const p=free.map(j=>h[j]/d);
  samples.push(p.every(v=>Number.isFinite(v)&&Math.abs(v)<=e)?p:null);
 }
 return samples;
}
function drawBoundaryLines(ctx,chart,rect,scale,view){
 for(const line of LINES){
  const selected=state.selected===line.label;ctx.beginPath();let drawing=false;
  for(const p of boundaryLineInChart(line,chart)){
   if(!p){drawing=false;continue}
   const q=rotate(p,view),x=rect.width/2+q[0]*scale,y=rect.height/2-q[1]*scale+8;
   if(drawing)ctx.lineTo(x,y);else{ctx.moveTo(x,y);drawing=true}
  }
  ctx.lineJoin="round";ctx.lineCap="round";
  ctx.strokeStyle="#082820";ctx.globalAlpha=.9;ctx.lineWidth=selected?6:4;ctx.stroke();
  ctx.strokeStyle=selected?"#f4e7dc":line.color;ctx.globalAlpha=selected?1:.9;ctx.lineWidth=selected?4:2.4;ctx.shadowColor=selected?"#f4e7dc":line.color;ctx.shadowBlur=selected?5:2;ctx.stroke();ctx.shadowBlur=0;
 }
 ctx.globalAlpha=1;
}
function drawChart(chart){
 const canvas=canvases[chart],ctx=canvas.getContext("2d",{alpha:false}),rect=canvas.getBoundingClientRect(),d=Math.min(devicePixelRatio||1,2),view=state.views[chart];canvas.width=Math.round(rect.width*d);canvas.height=Math.round(rect.height*d);ctx.setTransform(d,0,0,d,0,0);ctx.fillStyle="#123c35";ctx.fillRect(0,0,rect.width,rect.height);
 ctx.strokeStyle="#b6dbca14";for(let x=0;x<rect.width;x+=35){ctx.beginPath();ctx.moveTo(x,0);ctx.lineTo(x,rect.height);ctx.stroke()}for(let y=0;y<rect.height;y+=35){ctx.beginPath();ctx.moveTo(0,y);ctx.lineTo(rect.width,y);ctx.stroke()}
 const scale=Math.min(rect.width,rect.height)*.16*view.zoom,light=rotate([-.4,.65,.8],view),points=state.points[chart].map(({p,n})=>{const q=rotate(p,view),nr=rotate(n,view),shade=Math.max(0,nr[0]*light[0]+nr[1]*light[1]+nr[2]*light[2]);return{x:rect.width/2+q[0]*scale,y:rect.height/2-q[1]*scale+8,z:q[2],shade}}).sort((a,b)=>a.z-b.z);
 const layer=document.createElement("canvas"),lc=layer.getContext("2d");layer.width=canvas.width;layer.height=canvas.height;lc.setTransform(d,0,0,d,0,0);
 for(const p of points){lc.globalAlpha=.34+.42*p.shade;lc.fillStyle="#deef68";const s=4.0+1.6*p.shade;lc.fillRect(p.x-s/2,p.y-s/2,s,s)}
 ctx.save();ctx.globalAlpha=.72;ctx.filter="blur(.65px)";ctx.drawImage(layer,0,0,rect.width,rect.height);ctx.globalAlpha=.58;ctx.filter="none";ctx.drawImage(layer,0,0,rect.width,rect.height);ctx.restore()
 drawBoundaryLines(ctx,chart,rect,scale,view);
}
function drawAll(){canvases.forEach((_,i)=>drawChart(i))}
function planeSvg(family,l,r){
 const ys=[36,72,108],lx=105,rx=330;let lines="",dots="";
 for(let i=1;i<=3;i++)for(let j=1;j<=3;j++){const info=endpoints.get([`${l}${i}`,`${r}${j}`].sort().join("-"));lines+=`<line class="grid-line" data-line="${info.label}" x1="${lx}" y1="${ys[i-1]}" x2="${rx}" y2="${ys[j-1]}" stroke=""/>`}
 for(let i=1;i<=3;i++){dots+=`<circle class="dot" cx="${lx}" cy="${ys[i-1]}" r="3"/><text class="dotlabel" x="${lx-12}" y="${ys[i-1]+3}" text-anchor="end">${l}${i}</text><circle class="dot" cx="${rx}" cy="${ys[i-1]}" r="3"/><text class="dotlabel" x="${rx+12}" y="${ys[i-1]+3}">${r}${i}</text>`}
 return `<svg viewBox="0 0 410 145">${lines}${dots}</svg>`;
}
function build(){
 document.querySelector("#line-atlas").innerHTML=LINES.map(x=>"<button class=\"chip\" data-line=\""+x.label+"\" style=\"--line-color:"+x.color+"\"><b>"+x.label[0]+"<sub>"+x.label[1]+"</sub></b><span>"+x.from+"—"+x.to+"</span></button>").join("");
 document.querySelectorAll("[data-line]").forEach(el=>{el.addEventListener("mouseenter",()=>select(el.dataset.line));el.addEventListener("click",()=>{state.selected=el.dataset.line;select(el.dataset.line)})});
}
function select(label){
 const x=LINES.find(v=>v.label===label);if(!x)return;
 document.querySelectorAll("[data-line]").forEach(el=>el.classList.toggle("selected",el.dataset.line===label));
 drawAll();
}
function clear(){state.selected=null;document.querySelectorAll(".selected").forEach(x=>x.classList.remove("selected"));const d=document.querySelector("#line-detail");d.querySelector("small").textContent="Hover over a line";d.querySelector("strong").textContent="27 marked limits";d.querySelector("p").textContent="Each family contributes nine lines to one coordinate plane.";drawAll()}
function schedule(){updateEquation();document.querySelector("#surface-status").textContent="updating…";if(!timer)timer=setTimeout(sample,45)}
function updateRhoUI(){rhoInput.value=state.rho;document.querySelector("#rho-output").textContent=rhoZero()?"0 (limit)":state.rho.toFixed(3)}
function setRho(value){state.rho=Math.max(0,Math.min(1,value));updateRhoUI();document.querySelectorAll("[data-rho]").forEach(b=>b.classList.toggle("active",+b.dataset.rho===state.rho));schedule()}
const playing=new Map();
const playLabel=key=>({rho:"ρ",lambda:"λ",mu:"μ",nu:"ν"}[key]);
function stopPlay(key){const item=playing.get(key);if(!item)return;cancelAnimationFrame(item.frame);playing.delete(key);const b=document.querySelector("[data-play="+key+"]");b.textContent="▶ Play "+playLabel(key);b.classList.remove("playing");schedule()}
function togglePlay(key){
 if(playing.has(key)){stopPlay(key);return}
 const button=document.querySelector("[data-play="+key+"]");button.textContent="❚❚ Pause "+playLabel(key);button.classList.add("playing");
 const input=document.querySelector("#"+key),min=+input.min,max=+input.max;let direction=+input.value>=max-.001?-1:1,last=performance.now();
 const tick=now=>{let value=+input.value+direction*(max-min)*(now-last)/7000;last=now;if(value>=max){value=max;direction=-1}else if(value<=min){value=min;direction=1}input.value=value;input.dispatchEvent(new Event("input"));const item=playing.get(key);if(item)item.frame=requestAnimationFrame(tick)};playing.set(key,{frame:requestAnimationFrame(tick)});
}
document.querySelectorAll("[data-play]").forEach(b=>b.addEventListener("click",()=>togglePlay(b.dataset.play)));
["lambda","mu","nu"].forEach(key=>document.querySelector("#"+key).addEventListener("pointerdown",()=>stopPlay(key)));
rhoInput.addEventListener("pointerdown",()=>stopPlay("rho"));
rhoInput.addEventListener("input",e=>setRho(+e.target.value));
document.querySelectorAll("[data-rho]").forEach(b=>b.addEventListener("click",()=>setRho(+b.dataset.rho)));
["lambda","mu","nu"].forEach(k=>document.querySelector("#"+k).addEventListener("input",e=>{state[k]=+e.target.value;document.querySelector("#"+k+"-output").textContent=state[k].toFixed(2);schedule()}));

canvases.forEach((canvas,chart)=>{let dragging=false,last=[0,0];canvas.addEventListener("pointerdown",e=>{dragging=true;last=[e.clientX,e.clientY];canvas.setPointerCapture(e.pointerId)});canvas.addEventListener("pointermove",e=>{if(!dragging)return;const v=state.views[chart];v.yaw+=(e.clientX-last[0])*.008;v.pitch=Math.max(-1.35,Math.min(1.35,v.pitch+(e.clientY-last[1])*.008));last=[e.clientX,e.clientY];drawChart(chart)});canvas.addEventListener("pointerup",()=>dragging=false);canvas.addEventListener("wheel",e=>{e.preventDefault();const v=state.views[chart];v.zoom=Math.max(.55,Math.min(2.1,v.zoom*Math.exp(-e.deltaY*.001)));drawChart(chart)},{passive:false})});
addEventListener("resize",drawAll);
build();updateRhoUI();updateEquation();sample();
