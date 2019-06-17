import json
import io
import sys # for initializing min value by the max value
from collections import deque


class Graph:

    def __init__(self, fileName):
        self.myNodesToIndex = {}
        self.myNodesToId = {}
        self.myLinksRawData = []
        self.myNodesRawData = []
        self.numNode = 0
        self.numLinks = 0
        self.maxLevel = 0
        self.myMST = []
        self.marked = []
        self.myAdjMST = {}
        self.iterateData = []

        with open(fileName) as jsonData:
            rawData = json.load(jsonData)
            for node in rawData["nodes"]:
                # print(str(node["id"]) + " " + str(self.numNode) + "\n")
                self.myNodesToId[self.numNode] = node["label"]
                self.myNodesToIndex[node["label"]] = self.numNode
                self.myNodesRawData.append({
                    "id": node["label"],
                    "level": -1,
                    "isRootNode": False
                })
                self.numNode += 1

            self.myLinks = [[0 for column in range(self.numNode)] for row in range(self.numNode)]

            for link in rawData["links"]:
                # print(str(link) + "\n")
                # myLinks[link["source"]].append([link["target"], link["value"]])
                self.myLinks[self.myNodesToIndex[link["source"]]][self.myNodesToIndex[link["target"]]] = link["value"]
                self.myLinks[self.myNodesToIndex[link["target"]]][self.myNodesToIndex[link["source"]]] = link["value"]
                link["isTreeLink"] = False
                self.myLinksRawData.append(link)
                self.numLinks += 1

    def printMST(self):
        print("Edge \tWeight")
        for i in range(1, self.numNode):
            print(self.myMST[i], "-", i, "\t", self.myLinks[i][self.myMST[i]])

    def minKey(self, key, mstSet):
        min = sys.maxsize
        for vertexIndex in range(0, self.numNode):
            if key[vertexIndex] < min and mstSet[vertexIndex] == False:
                min = key[vertexIndex]
                min_index = vertexIndex

        return min_index

    def primMST(self):

        # Key values used to pick minimum weight edge in cut
        key = [sys.maxsize] * self.numNode
        self.myMST = [None] * self.numNode  # Array to store constructed MST
        # Make key 0 so that this vertex is picked as first vertex
        key[0] = 0
        mstSet = [False] * self.numNode

        self.myMST[0] = -1  # First node is always the root of

        for cout in range(self.numNode):

            # Pick the minimum distance vertex from
            # the set of vertices not yet processed.
            # u is always equal to src in first iteration
            u = self.minKey(key, mstSet)

            # Put the minimum distance vertex in
            # the shortest path tree
            mstSet[u] = True

            # Update dist value of the adjacent vertices
            # of the picked vertex only if the current
            # distance is greater than new distance and
            # the vertex in not in the shotest path tree
            for v in range(self.numNode):
                # graph[u][v] is non zero only for adjacent vertices of m
                # mstSet[v] is false for vertices not yet included in MST
                # Update the key only if graph[u][v] is smaller than key[v]
                if self.myLinks[u][v] > 0 and mstSet[v] == False and key[v] > self.myLinks[u][v]:
                    key[v] = self.myLinks[u][v]
                    self.myMST[v] = u

        # print("MST Generated: " + str(self.myMST))
        # print("MST Assist: ")
        # self.printMST()

        # self.myAdjMST = [[0 for i in range(self.numNode)] for j in range(0)]
        for i in range(0, self.numNode):
            self.myAdjMST[i] = []

        for i in range(1, self.numNode):
            self.myAdjMST[i].append(self.myMST[i]) #also inversed
            self.myAdjMST[self.myMST[i]].append(i)

        dummyLinks = []
        for i in range(1, self.numNode):
            # print(self.myMST[i], "-", i, "\t", self.myLinks[i][self.myMST[i]])
            # TODO: optimize this part
            dummyLinks.append({
                "source": self.myNodesToId[self.myMST[i]],
                "target": self.myNodesToId[i],
                "value": self.myLinks[i][self.myMST[i]]
            })
            # print(dummyLinks)

        # print("TreeLinks: ")
        # print(dummyLinks)

        for link in self.myLinksRawData:
            isFoundTreeLink = False
            for treeLink in dummyLinks:
                if (link["source"] == treeLink["source"] and link["target"] == treeLink["target"]):
                    link["isTreeLink"] = True
                    isFoundTreeLink = True
            if (not isFoundTreeLink):
                link["isTreeLink"] = False

    def markLevelInfo(self):
        # print("Given MST Adjmat: " + str(self.myAdjMST))
        level = [-1] * self.numNode
        self.marked = [False] * self.numNode

        que = deque([])

        #root
        que.append(0)

        level[0] = 0
        self.marked[0] = True

        while(len(que) != 0) :
            node = que.popleft()
            for i in range(len(self.myAdjMST[node])):
                neighbor = self.myAdjMST[node][i]

                if (not self.marked[neighbor]):
                    que.append(neighbor)
                    level[neighbor] = level[node] + 1
                    self.marked[neighbor] = True
                    if level[neighbor] > self.maxLevel:
                        self.maxLevel = level[neighbor]

        for node in self.myNodesRawData:
            # print("level: " + str(level[self.myNodesToIndex[node["id"]]]))
            node["level"] = level[self.myNodesToIndex[node["id"]]]

    def generateGroupingData(self):
        groupNodes = []
        dummyGroupNodes = []
        groupLinks = []
        groupNodesHashSet = {}
        dummyGroupNodesHashSet = {}

        for i in range(0, self.maxLevel):
            groupNodes.clear()
            groupLinks.clear()
            dummyGroupNodes.clear()
            groupNodesHashSet.clear()
            dummyGroupNodesHashSet.clear()

            # print("clear: ")
            # print(dummyGroupNodesHashSet)
            for node in self.myNodesRawData:
                if node["level"] == i or node["level"] == i + 1:
                    # print(i)
                    # print(node["level"])
                    dummyGroupNodes.append(node)
                    # groupNodes.append(node)
                    groupNodesHashSet[node["id"]] = True

            print("GroupNodes: ")
            print(groupNodesHashSet)

            outputLinks = []
            for i in range(1, self.numNode):
                # print(self.myMST[i], "-", i, "\t", self.myLinks[i][self.myMST[i]])
                outputLinks.append({
                    "source": self.myNodesToId[i],  # this is reversed
                    "target": self.myNodesToId[self.myMST[i]],
                    "value": self.myLinks[i][self.myMST[i]],
                })

            for link in outputLinks: # self.myLinksRawData:
                # if (self.myNodesRawData[self.myNodesToIndex[link["source"]]]["level"] == i or self.myNodesRawData[self.myNodesToIndex[link["source"]]]["level"] == i + 1) and (self.myNodesRawData[self.myNodesToIndex[link["target"]]]["level"] == i or self.myNodesRawData[self.myNodesToIndex[link["target"]]]["level"] == i + 1) and link["isTreeLink"]:
                if (link["source"] in groupNodesHashSet) and (link["target"] in groupNodesHashSet):
                    print("Found in set:")
                    print(link)
                    # if link["isTreeLink"]:
                        # print(link)
                    groupLinks.append(link)
                    link["isTreeLink"] = True
                    dummyGroupNodesHashSet[link["source"]] = True
                    dummyGroupNodesHashSet[link["target"]] = True
            # this block have some bug TODO: bug fix
            # TODO: bug fix! unidentified nodes links issue cause links all blank
            # for node in dummyGroupNodes:
            #    for link in groupLinks:
            #        if node["id"] == link["source"]["id"] or node["id"] == link["target"]["id"]:
            #            groupNodes.append(node)

            print("DummyNodeHashSet: ")
            print(dummyGroupNodesHashSet)

            for key in dummyGroupNodesHashSet:
                groupNodes.append(self.myNodesRawData[self.myNodesToIndex[key]])

            self.iterateData.append({
                "nodes": groupNodes.copy(),
                "links": groupLinks.copy()
            })




    def dumpJSON(self, onlyMSTLinks = False):
        print("Dumping...")
        self.markLevelInfo()
        outputLinks = []
        if onlyMSTLinks:
            for i in range(1, self.numNode):
                # print(self.myMST[i], "-", i, "\t", self.myLinks[i][self.myMST[i]])
                outputLinks.append({
                    "source": self.myNodesToId[self.myMST[i]],
                    "target": self.myNodesToId[i],
                    "value": self.myLinks[i][self.myMST[i]],
                    "isTreeLink": True
                })
        else:
            dummyLinks = []
            for i in range(1, self.numNode):
                # print(self.myMST[i], "-", i, "\t", self.myLinks[i][self.myMST[i]])
                # TODO: optimize this part
                dummyLinks.append({
                    "source": self.myNodesToId[self.myMST[i]],
                    "target": self.myNodesToId[i],
                    "value": self.myLinks[i][self.myMST[i]]
                })
                #print(dummyLinks)

            # print("TreeLinks: ")
            # print(dummyLinks)

            for link in self.myLinksRawData:
                isFoundTreeLink = False
                for treeLink in dummyLinks:
                    if (link["source"] == treeLink["source"] and link["target"] == treeLink["target"]):
                        outputLinks.append({
                            "source": link["source"],
                            "target": link["target"],
                            "value": link["value"],
                            "isTreeLink": True
                        })
                        isFoundTreeLink = True
                if (not isFoundTreeLink) :
                    outputLinks.append({
                        "source": link["source"],
                        "target": link["target"],
                        "value": link["value"],
                        "isTreeLink": False
                    })

        outputNodes = []

        for node in self.myNodesRawData:
            if self.myNodesToIndex[node["id"]] == 0:
                outputNodes.append({
                    "id": node["id"],
                    "level": node["level"],
                    "isRootNode": True
                })
            else:
                outputNodes.append({
                    "id": node["id"],
                    "level": node["level"],
                    "isRootNode": False
                })

        '''
        try:
            to_unicode = unicode
        except NameError:
            to_unicode = str

        with io.open('data.json', 'w', encoding='utf8') as outfile:
            str_ = json.dumps(data, indent=4, sort_keys=False, separators=(',', ': '), ensure_ascii=False)
            outfile.write(str_)
        '''

        # import data
        data = {
            "nodes": outputNodes,
            "links": outputLinks,
            'metadata': {
                "NumNodes": self.numNode,
                "NumLinks": self.numLinks,
                "MaxLevel": self.maxLevel
            }
        }

        with io.open('processedGraph.json', 'w', encoding='utf8') as outfile:
            str_ = json.dumps(data, indent=4, sort_keys=False, separators=(',', ': '), ensure_ascii=False)
            outfile.write(str_)

        print("Finish Dumping!")

    def dumpJSONWithSep(self):
        print("Dumping...")
        self.markLevelInfo()
        self.generateGroupingData()

        outputLinks = []
        for i in range(1, self.numNode):
            # print(self.myMST[i], "-", i, "\t", self.myLinks[i][self.myMST[i]])
            outputLinks.append({
                "source": self.myNodesToId[i], # this is reversed
                "target": self.myNodesToId[self.myMST[i]],
                "value": self.myLinks[i][self.myMST[i]],
            })

        outputNodes = []
        for node in self.myNodesRawData:
            if self.myNodesToIndex[node["id"]] == 0:
                outputNodes.append({
                    "id": node["id"],
                    "level": node["level"],
                    "isRootNode": True
                })
            else:
                outputNodes.append({
                    "id": node["id"],
                    "level": node["level"],
                    "isRootNode": False
                })

        # import data
        data = {
            "nodes": outputNodes,
            #"allLinks": self.myLinksRawData,
            "links": outputLinks,
            'metadata': {
                "NumNodes": self.numNode,
                "NumLinks": self.numLinks,
                "MaxLevel": self.maxLevel
            },
            #"iterateData": self.iterateData
        }

        with io.open('processedPower.json', 'w', encoding='utf8') as outfile:
            str_ = json.dumps(data, indent=4, sort_keys=False, separators=(',', ': '), ensure_ascii=False)
            outfile.write(str_)

        print("Finish Dumping!")


if __name__ == '__main__':
    fileName = "power.json"
    graph = Graph(fileName)
    graph.primMST()
    graph.dumpJSONWithSep()

