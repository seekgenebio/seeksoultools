import sys
from multiprocessing import Pipe, Process, Queue
import multiprocessing.connection
import io
import dnaio
from xopen import xopen

class Reader(Process):
    """
    读取paired fastq
    """
    def __init__(self, file1:str, file2:str, connections:multiprocessing.connection, queue:Queue, buffer_size:int):
        """初始化读入进程

        Args:
            file1: read1 fastq file
            file2: read2 fastq file
            connections: 
            queue: 存放空闲工作进程index的队列
            buffer_size: 读取fastq时的buffer大小
        
        Returns:
            Reader对象

        """
        super().__init__()
        self.file1 = file1
        self.file2 = file2
        self.connections = connections
        self.queue = queue
        self.buffer_size = buffer_size

    def run(self):
        try:
            chunk_index = 0
            for file1, file2 in zip(self.file1, self.file2):
                with xopen(file1, 'rb') as f1:
                    with xopen(file2, 'rb') as f2:
                        for (chunk1, chunk2) in dnaio.read_paired_chunks(f1, f2, self.buffer_size):
                            worker_index = self.queue.get()
                            pipe = self.connections[worker_index]
                            pipe.send(chunk_index)
                            pipe.send_bytes(chunk1)
                            pipe.send_bytes(chunk2)
                            chunk_index += 1
            for _ in range(len(self.connections)):
                worker_index = self.queue.get()
                self.connections[worker_index].send(-1)
        except Exception as e:
            for worker_index in range(len(self.connections)):
                self.connections[worker_index].send(-2)
            raise e

class Writer:
    """
    处理输出内容
    """
    def __init__(self, file1:str, file_multi:str, file2:str=None):
        """输出处理好的序列
        Args:
            file: read2 fastq file1
            file_muti: 不能确定barcode的fastq文件名称
        """
        self._file1 = file1
        self._fh1 = xopen(self._file1, mode='wb')

        self._file_multi = file_multi
        self._fh_multi = xopen(self._file_multi, mode='wb')

        self._file2 = file2
        if self._file2:
            self._fh2 = xopen(self._file2, mode='wb')
        self._chunks = dict()
        self._current_index = 0

    def write(self, data, index):
        self._chunks[index] = data
        while self._current_index in self._chunks:
            self._fh1.write(self._chunks[self._current_index][0])
            self._fh_multi.write(self._chunks[self._current_index][1])
            if self._file2:
                self._fh2.write(self._chunks[self._current_index][2])
            del self._chunks[self._current_index]
            self._current_index += 1

    def wrote_everything(self):
        return not self._chunks

    def close(self):
        self._fh1.close()
        self._fh_multi.close()
        if self._file2:
            self._fh2.close()

class Worker(Process):
    """工作进程类
    """
    def __init__(self, id_, read_pipe, write_pipe, need_work_queue, func, paired_out=False):
        """
        """
        super().__init__()
        self._id = id_
        self.read_pipe = read_pipe
        self.write_pipe = write_pipe
        self.need_work_queue = need_work_queue
        self.func = func
        self.paired_out = paired_out

    def run(self):
        try:
            while True:
                self.need_work_queue.put(self._id)
                chunk_index = self.read_pipe.recv()
                if chunk_index == -1:
                    break
                elif chunk_index == -2:
                    e, tb_str = self.read_pipe.recv()
                    raise e
                data = self.read_pipe.recv_bytes()
                input = io.BytesIO(data)
                data = self.read_pipe.recv_bytes()
                input2 = io.BytesIO(data)
                tmp = io.BytesIO()
                tmp_multi = io.BytesIO()
                if self.paired_out:
                    tmp2 = io.BytesIO()
                    _ = self.func(fq1=input, fq2=input2, fq_out=tmp, fqout_multi=tmp_multi, fq_out2=tmp2)
                else:
                    _ = self.func(fq1=input, fq2=input2, fq_out=tmp, fqout_multi=tmp_multi)
                self.write_pipe.send(chunk_index)
                self.write_pipe.send_bytes(tmp.getvalue())
                self.write_pipe.send_bytes(tmp_multi.getvalue())
                if self.paired_out:
                    self.write_pipe.send_bytes(tmp2.getvalue())
                self.write_pipe.send(_)
            self.write_pipe.send(-1)
        except Exception as e:
            self.write_pipe.send(-2)
            raise e


class Pipeline:
    def __init__(self, func, fq1, fq2, fqout1, fqout_multi, core, stat=None, fqout2=None, buffer_size=16*1024**2):
        self.n_workers = core
        self.fq1 = fq1
        self.fq2 = fq2
        self.fqout1 = fqout1
        self.fqout_multi = fqout_multi
        self.fqout2 = fqout2
        self.buffer_size = buffer_size
        self.need_work_queue = Queue()
        self.func = func
        self.paired_out = False
        if self.fqout2:
            self.paired_out = True
        self.stat = stat

    def run(self):
        # start reader process
        reader_connections = [Pipe(duplex=False) for _ in range(self.n_workers)]
        _pipes, _conn = zip(*reader_connections)
        _reader_process = Reader(self.fq1, self.fq2, _conn, self.need_work_queue, self.buffer_size)
        _reader_process.daemon = True
        _reader_process.start()

        # start worker processes
        self.workers = []
        self.connections = []
        self.writer = Writer(self.fqout1, self.fqout_multi, self.fqout2)
        for index in range(self.n_workers):
            conn_r, conn_w = Pipe(duplex=False)
            self.connections.append(conn_r)
            worker = Worker(index, _pipes[index], conn_w, self.need_work_queue,
                            self.func, self.paired_out)
            worker.daemon = True
            worker.start()
            self.workers.append(worker)

        # write output
        while self.connections:
            ready_connections = multiprocessing.connection.wait(self.connections)
            for connection in ready_connections:
                chunk_index = connection.recv()
                if chunk_index == -1:
                    self.connections.remove(connection)
                    continue
                elif chunk_index == -2:
                    sys.stderr.write('err!!!\n')
                
                if self.paired_out:
                    data1 = connection.recv_bytes()
                    data_multi = connection.recv_bytes()
                    data2 = connection.recv_bytes()
                    self.writer.write([data1, data_multi, data2], chunk_index)
                else:
                    data1 = connection.recv_bytes()
                    data_multi = connection.recv_bytes()
                    self.writer.write([data1, data_multi], chunk_index)
                _stat = connection.recv()
                self.stat.update(**_stat)
        assert self.writer.wrote_everything()
        for w in self.workers:
            w.join()
        _reader_process.join()
        self.writer.close()
