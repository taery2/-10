"""ЛР 10, вариант 4. Кратчайшие пути, Дейкстра и критический путь."""
import random, time, heapq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
INF = float("inf")

def make_graph(n, m, seed=42):
    rng = random.Random(seed)
    edges = set()
    while len(edges) < m:
        u = rng.randrange(n); v = rng.randrange(n)
        if u != v: edges.add((min(u,v), max(u,v)))
    adj = [[] for _ in range(n)]
    for u,v in edges:
        w = rng.randint(1,100)
        adj[u].append((v,w)); adj[v].append((u,w))
    return adj

def dijkstra_heap(adj, source):
    dist = [INF]*len(adj); parent = [-1]*len(adj)
    dist[source] = 0; heap = [(0,source)]
    while heap:
        d,v = heapq.heappop(heap)
        if d != dist[v]: continue
        for u,w in adj[v]:
            nd = d+w
            if nd < dist[u]:
                dist[u] = nd; parent[u] = v
                heapq.heappush(heap,(nd,u))
    return dist,parent

def dijkstra_array(adj, source):
    n = len(adj); dist = [INF]*n; used = [False]*n; parent = [-1]*n
    dist[source] = 0
    for _ in range(n):
        v = min((i for i in range(n) if not used[i]), key=lambda i:dist[i], default=-1)
        if v == -1 or dist[v] == INF: break
        used[v] = True
        for u,w in adj[v]:
            if dist[v]+w < dist[u]:
                dist[u] = dist[v]+w; parent[u] = v
    return dist,parent

def restore_path(parent, target):
    path=[]; v=target
    while v != -1:
        path.append(v); v=parent[v]
    return list(reversed(path))

def topo(tasks):
    indeg={v:len(pred) for v,(dur,pred) in tasks.items()}
    succ={v:[] for v in tasks}
    for v,(dur,pred) in tasks.items():
        for p in pred: succ[p].append(v)
    q=[v for v,d in indeg.items() if d==0]; q.sort()
    order=[]
    while q:
        v=q.pop(0); order.append(v)
        for u in succ[v]:
            indeg[u]-=1
            if indeg[u]==0: q.append(u)
    if len(order)!=len(tasks): raise ValueError("Цикл в сетевом графике")
    return order,succ

def critical_path(tasks):
    order,succ=topo(tasks); ES={}; EF={}
    for v in order:
        dur,pred=tasks[v]; ES[v]=max((EF[p] for p in pred),default=0); EF[v]=ES[v]+dur
    T=max(EF.values()); LS={}; LF={}
    for v in reversed(order):
        dur=tasks[v][0]; LF[v]=min((LS[u] for u in succ[v]),default=T); LS[v]=LF[v]-dur
    R={v:LS[v]-ES[v] for v in tasks}
    return T,ES,EF,LS,LF,R,[v for v in order if R[v]==0]

if __name__ == "__main__":
    print("ЛР 10. Булова Екатерина, БИ 2-1, вариант 4")
    print("n       m       Дейкстра-массив, с   Дейкстра-куча, с")
    ns=[]; ta=[]; th=[]
    for n in (1000,2000,4000):
        m=5*n
        adj=make_graph(n,m,seed=200+n)
        t=time.perf_counter(); da,_=dijkstra_array(adj,0); a=time.perf_counter()-t
        t=time.perf_counter(); dh,ph=dijkstra_heap(adj,0); h=time.perf_counter()-t
        assert da==dh
        print(f"{n:<8}{m:<8}{a:<22.6f}{h:.6f}")
        ns.append(n); ta.append(a); th.append(h)

    # 11 работ проекта (A–K), зависимости образуют DAG.
    tasks={
      "A":(3,[]), "B":(5,["A"]), "C":(4,["A"]), "D":(6,["B"]),
      "E":(7,["B","C"]), "F":(3,["C"]), "G":(5,["D","E"]),
      "H":(4,["E","F"]), "I":(6,["G"]), "J":(2,["H","I"]),
      "K":(3,["J"])
    }
    T,ES,EF,LS,LF,R,crit=critical_path(tasks)
    print(f"\\nКритический путь: {' -> '.join(crit)}; длительность проекта: {T}")
    print("Работа  t  ES  EF  LS  LF  Резерв")
    for v in tasks:
        print(f"{v:>4} {tasks[v][0]:>3} {ES[v]:>3} {EF[v]:>3} {LS[v]:>3} {LF[v]:>3} {R[v]:>5}")

    # Доп. К4: количество кратчайших путей.
    # Малый граф с двумя равными по длине маршрутами 0-1-3 и 0-2-3.
    test=[[(1,2),(2,2)],[(0,2),(3,3)],[(0,2),(3,3)],[(1,3),(2,3)]]
    dist=[INF]*4; count=[0]*4; dist[0]=0; count[0]=1; heap=[(0,0)]
    while heap:
        d,v=heapq.heappop(heap)
        if d!=dist[v]: continue
        for u,w in test[v]:
            nd=d+w
            if nd<dist[u]:
                dist[u]=nd; count[u]=count[v]; heapq.heappush(heap,(nd,u))
            elif nd==dist[u]:
                count[u]+=count[v]
    assert dist[3]==5 and count[3]==2
    print(f"\\nК4. Кратчайшее расстояние 0->3: {dist[3]}, число кратчайших путей: {count[3]}")

    plt.figure()
    plt.plot(ns,ta,marker="o",label="Поиск минимума в массиве")
    plt.plot(ns,th,marker="s",label="Приоритетная очередь")
    plt.xlabel("Число вершин n"); plt.ylabel("Время, с")
    plt.title("ЛР 10. Сравнение реализаций Дейкстры"); plt.grid(True); plt.legend()
    plt.tight_layout(); plt.savefig("results/dijkstra_time.png",dpi=160)
