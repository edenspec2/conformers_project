from contextlib import contextmanager
import threading
import _thread

class TimeLimiterForLoops():
    def __init__(self):
        self.skipped_indices=[]

    @contextmanager
    def time_limit(self, seconds, loop_index=0):
        timer=threading.Timer(seconds, lambda: _thread.interrupt_main())
        timer.start()
        try:
            yield
        except KeyboardInterrupt:
            print('Skipped Entry')
            self.skipped_indices.append(loop_index)
        finally:
            timer.cancel()

@contextmanager
def time_limit(seconds):
    timer=threading.Timer(seconds, lambda: _thread.interrupt_main())
    timer.start()
    try:
        yield
    except KeyboardInterrupt:
        print('skipped run')
    finally:
        # if the action ends in specified time, timer is canceled
        timer.cancel()


if __name__=='__main__':
    import time
    with time_limit(ENTRY_PROCESSING_TIME_LIMIT):
        for i in range(20):
            print(i)
            time.sleep(1)    
    print('script loaded successfully')


